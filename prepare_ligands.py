# Имя файла: prepare_ligands_flexible.py
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import os
import sys

def prepare_ligands_for_interformer(input_csv_path, output_sdf_path, smiles_column='ligand', max_mol_weight=500):
    """
    Гибкая функция для подготовки лигандов.
    Принимает пути в качестве аргументов.
    """
    print("--- Запуск гибкого скрипта подготовки лигандов ---")
    
    output_dir = os.path.dirname(output_sdf_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Создана директория: {output_dir}")

    try:
        df = pd.read_csv(input_csv_path, delimiter=',')
        if smiles_column not in df.columns:
            print(f"КРИТИЧЕСКАЯ ОШИБКА: Колонка '{smiles_column}' не найдена в файле {input_csv_path}.")
            print(f"Найденные колонки: {list(df.columns)}")
            return
            
    except FileNotFoundError:
        print(f"КРИТИЧЕСКАЯ ОШИБКА: Файл не найден: {input_csv_path}")
        return
    except Exception as e:
        print(f"Произошла ошибка при чтении CSV: {e}")
        return

    writer = Chem.SDWriter(output_sdf_path)
    print(f"Начинаем обработку {len(df)} лигандов из колонки '{smiles_column}'...")
    
    success_count, filtered_out_count, error_count = 0, 0, 0

    for index, row in df.iterrows():
        smiles = row.get(smiles_column)
        if pd.isna(smiles):
            error_count += 1
            continue
        
        mol = Chem.MolFromSmiles(str(smiles).strip())
        
        if mol is None:
            error_count += 1
            continue

        if Descriptors.MolWt(mol) > max_mol_weight:
            filtered_out_count += 1
            continue
            
        try:
            mol_h = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            AllChem.EmbedMolecule(mol_h, params)
            AllChem.UFFOptimizeMolecule(mol_h)
            mol_h.SetProp("_Name", f"ligand_source_row_{index+2}")
            writer.write(mol_h)
            success_count += 1
        except (RuntimeError, ValueError):
            error_count += 1
            continue

    writer.close()
    
    print("\n--- Отчет ---")
    print(f"Успешно сохранено: {success_count} молекул")
    print(f"Отфильтровано по массе (> {max_mol_weight} Da): {filtered_out_count} молекул")
    print(f"Пропущено из-за ошибок: {error_count} молекул")
    print(f"Результат сохранен в файл: {output_sdf_path}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Ошибка: неверное количество аргументов.")
        print("Пример использования: python prepare_ligands_flexible.py <путь_к_входному_csv> <путь_к_выходному_sdf>")
        sys.exit(1)
    
    input_csv = sys.argv[1]
    output_sdf = sys.argv[2]
    prepare_ligands_for_interformer(input_csv, output_sdf)
