#!/bin/bash

#=========================================================================================
# SBATCH –- это инструкции для планировщика Slurm
#=========================================================================================
#SBATCH --job-name=interformer_docking
#SBATCH --partition=aichem
#SBATCH --gres=gpu:1
#SBATCH --mem=50G
#SBATCH --time=16:00:00
#SBATCH --output=slurm_output_%j.log
#SBATCH --error=slurm_errors_%j.err

#=========================================================================================
# КОНФИГУРАЦИЯ ЗАПУСКА (редактируйте здесь)
#=========================================================================================

# --- Основные параметры ---
PDB_ID="1tqn"
WORK_PATH="vs"
DOCKING_PATH="energy_VS"
LIGAND_FILE="CYP2d6_ki_df"
OMP_THREADS="1,1" # Для параллельного режима лучше ставить 1

# --- Параметры производительности ---
# Количество параллельных процессов для докинга.
# Это значение будет передано и в Slurm (--cpus-per-task), и в Python-скрипт.
NUM_WORKERS=10
#SBATCH --cpus-per-task=${NUM_WORKERS}

#=========================================================================================
# Настройка окружения
#=========================================================================================
echo "--- Starting job $SLURM_JOB_ID on host $SLURMD_NODENAME ---"
echo "Requesting ${NUM_WORKERS} CPUs for parallel processing."

echo "Activating conda environment..."
eval "$(conda shell.bash hook)"
conda activate interformer_env # <-- Убедитесь, что здесь имя вашего окружения

if [ $? -ne 0 ]; then
    echo "Failed to activate conda environment."
    exit 1
fi
echo "Conda environment activated successfully."

#=========================================================================================
# Запуск основного Python-скрипта
#=========================================================================================
echo "--- Starting master_pipeline.py in parallel mode ---"
python master_pipeline.py \
    -p ${PDB_ID} \
    -w ${WORK_PATH} \
    -d ${DOCKING_PATH} \
    -l ${LIGAND_FILE} \
    -o ${OMP_THREADS} \
    --parallel_workers ${NUM_WORKERS}

echo "Job finished with exit code $?"
