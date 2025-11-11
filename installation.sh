git clone https://github.com/tencent-ailab/Interformer.git
cd Interformer

conda env create --file environment.yml

if you have this problem:
OSError: [Errno 122] Disk quota exceeded: '/nfs/home/pshestun/miniconda3/pkgs/cache/ffeee55f.b418.tmp'

or mistake like this
use: 
conda clean --all

if it can't help, use then:

conda env list
conda remove --name <<YOUR_ENV_NAME>> --all


conda activate interformer

conda install -c conda-forge plip 
! We have any problem with original methods from github, and use this alternative

after all we need to use:
cd docking && python setup.py install && cd ../


conda install -c conda-forge openbabel

! And now we need to install checkpoints of model:

wget https://zenodo.org/records/15694429/files/checkpoints.zip?download=1
unzip checkpoints.zip

and make new file "run_interformer.sh":
nano run_interformer.sh

#!/bin/bash

#==============================================================================
# Slurm-скрипт для полного запуска тестового пайплайна Interformer
# Автор: AI Assistant (на основе отладки с пользователем)
# Версия: 2.0 (Окончательно исправленная)
#==============================================================================

# --- Параметры для планировщика Slurm ---
#SBATCH --job-name=interformer_pipeline  # Имя задачи в очереди
#SBATCH --output=interformer_pipeline_%j.out # Файл для стандартного вывода
#SBATCH --error=interformer_pipeline_%j.err  # Файл для вывода ошибок
#SBATCH --time=01:00:00                # Максимальное время выполнения (1 час)
#SBATCH --nodes=1                      # Запрашиваем 1 вычислительный узел
#SBATCH --ntasks-per-node=1            # 1 задача на узел
#SBATCH --cpus-per-task=8              # Количество ядер CPU для задачи
#SBATCH --mem=32G                      # Объем оперативной памяти
#SBATCH --gres=gpu:1                   # Запрос одного GPU
#SBATCH --partition=aichem             # Указание раздела (очереди) с GPU

# --- Подготовка окружения ---
echo "========================================================================"
echo "Starting Slurm job ${SLURM_JOB_ID} on host ${HOSTNAME}"
date
echo "========================================================================"
set -e
set +u
source "/nfs/home/pshestun/miniconda3/etc/profile.d/conda.sh"
conda activate interformer
set -u
echo "Conda environment 'interformer' activated successfully."
echo ""

# --- Основная часть скрипта ---

# --- ЭТАП 1: Подготовка данных ---
echo "--- [Этап 1/4] Начало подготовки данных ---"
echo "1.0 Создание необходимых директорий..."
mkdir -p examples/ligand examples/uff examples/pocket examples/raw/pocket

echo "1.1 Протонирование лиганда..."
obabel examples/raw/2qbr_ligand.sdf -p 7.4 -O examples/ligand/2qbr_docked.sdf

echo "1.2 Генерация начальной конформации лиганда..."
python tools/rdkit_ETKDG_3d_gen.py examples/ligand/ examples/uff

echo "1.3 Предобработка белка..."
reduce examples/raw/2qbr.pdb > examples/raw/pocket/2qbr_reduce.pdb

echo "1.4 Извлечение кармана связывания..."
python tools/extract_pocket_by_ligand.py examples/raw/pocket/ examples/ligand/ 1

# =============================================================================
# ФИНАЛЬНОЕ ИСПРАВЛЕНИЕ: Указываем правильный путь к выходному файлу
# =============================================================================
echo "1.5 Перемещение файла кармана из 'examples/raw/pocket/output'..."
mv examples/raw/pocket/output/2qbr_pocket.pdb examples/pocket
echo "1.6 Очистка временной папки 'examples/raw/pocket/output'..."
rm -rf examples/raw/pocket/output
# =============================================================================

echo "--- [Этап 1/4] Подготовка данных успешно завершена. ---"
echo ""

# --- ЭТАП 2: Предсказание энергетических функций ---
echo "--- [Этап 2/4] Начало предсказания энергетических функций ---"
DOCK_FOLDER="energy_output"
if [ ! -d "checkpoints" ]; then
    echo "КРИТИЧЕСКАЯ ОШИБКА: Директория 'checkpoints' не найдена!"
    exit 1
fi
PYTHONPATH=interformer/ python inference.py -test_csv examples/demo_dock.csv \
-work_path examples/ \
-ensemble checkpoints/vо.2_energy_model \
-batch_size 1 \
-posfix *val_loss* \
-energy_output_folder $DOCK_FOLDER \
-reload \
-debug
echo "--- [Этап 2/4] Предсказание энергий успешно завершено. ---"
echo ""

# --- ЭТАП 3: Докинг (генерация поз) и обработка результатов ---
echo "--- [Этап 3/4] Начало докинга и обработки результатов ---"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "3.1 Запуск генерации поз (докинг) с ${OMP_NUM_THREADS} потоками..."
python docking/reconstruct_ligands.py -y --cwd $DOCK_FOLDER -y --find_all find
echo "3.2 Создание сводного CSV-файла..."
python docking/reconstruct_ligands.py --cwd $DOCK_FOLDER --find_all stat
echo "3.3 Объединение CSV-файлов..."
python docking/merge_summary_input.py $DOCK_FOLDER/ligand_reconstructing/stat_concated.csv examples/demo_dock.csv
echo "--- [Этап 3/4] Докинг и обработка результатов успешно завершены. ---"
echo ""

# --- ЭТАП 4: Оценка поз (PoseScore) и аффинности ---
echo "--- [Этап 4/4] Начало оценки поз и аффинности ---"
echo "4.1 Копирование поз для оценки..."
mkdir -p examples/infer
cp -r $DOCK_FOLDER/ligand_reconstructing/*.sdf examples/infer
echo "4.2 Очистка кеша 'tmp_beta'..."
rm -rf tmp_beta
echo "4.3 Запуск предсказания PoseScore и аффинности..."
PYTHONPATH=interformer/ python inference.py -test_csv examples/demo_dock.round0.csv \
-work_path examples/ \
-ligand_folder infer/ \
-ensemble checkpoints/v0.2_affinity_model/model* \
-gpus 1 \
-batch_size 20 \
-posfix *val_loss* \
--pose_sel True
echo "--- [Этап 4/4] Оценка успешно завершена. ---"
echo ""

# --- Завершение работы ---
echo "========================================================================"
echo "Все этапы пайплайна Interformer успешно завершены!"
echo "Финальные результаты находятся в файле: result/demo_dock.round0_ensemble.csv"
date
echo "========================================================================"
