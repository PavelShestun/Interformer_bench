import argparse
import logging
import os
import shutil
import subprocess
import sys
import pandas as pd
from rdkit import Chem
from multiprocessing import Pool, cpu_count

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(processName)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("master_pipeline.log"),
        logging.StreamHandler(sys.stdout)
    ]
)

class PipelineStepFailed(Exception):
    """Custom exception for failed pipeline steps."""
    pass

def run_command(command, log_file_path=None):
    """Executes a command and logs its output."""
    logging.info(f"Executing command: {' '.join(command)}")
    try:
        # For parallel execution, it's better to capture output to files
        with open(log_file_path, 'a') as log_file:
            result = subprocess.run(command, check=True, stdout=log_file, stderr=log_file, text=True)
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {' '.join(command)}. See log for details: {log_file_path}")
        # The error is already logged to the file, so we just raise
        raise PipelineStepFailed(f"Command failed: {' '.join(command)}")

def prepare_directories(work_path):
    """Creates the necessary directories for the pipeline."""
    for subdir in ["data", "raw", "raw/pocket", "ligand", "pocket", "uff", "infer", "result"]:
        os.makedirs(os.path.join(work_path, subdir), exist_ok=True)
    os.makedirs("result", exist_ok=True)


def copy_input_files(work_path, pdb, df_lig):
    """Copies input files to the appropriate directories."""
    logging.info("Copying input files...")
    shutil.copy(os.path.join(work_path, "data", f"{pdb}_ligand.pdb"), os.path.join(work_path, "raw", f"{pdb}_ligand.pdb"))
    shutil.copy(os.path.join(work_path, "data", f"{df_lig}.csv"), os.path.join(work_path, "raw", f"{df_lig}.csv"))
    shutil.copy(os.path.join(work_path, "data", f"{pdb}_meeko.pdb"), os.path.join(work_path, "raw", "pocket", f"{pdb}_meeko.pdb"))

def prepare_protein_and_ligand(work_path, pdb):
    """Prepares the protein and reference ligand."""
    logging.info("Preparing protein and reference ligand...")
    log_path = os.path.join("result", "setup_commands.log")
    run_command(["obabel", os.path.join(work_path, "raw", f"{pdb}_ligand.pdb"), "-p", "7.4", "-O", os.path.join(work_path, "ligand", f"{pdb}_docked.sdf")], log_path)
    run_command(["python", "tools/extract_pocket_by_ligand.py", os.path.join(work_path, "raw", "pocket/"), os.path.join(work_path, "ligand/"), "0"], log_path)

    source_path = os.path.join(work_path, "raw", "pocket", "output", f"{pdb}_pocket.pdb")
    destination_path = os.path.join(work_path, "pocket", f"{pdb}_pocket.pdb")

    if not os.path.exists(source_path):
        raise PipelineStepFailed(f"Source file for pocket not found after extraction: {source_path}")

    if os.path.exists(destination_path):
        os.remove(destination_path)
    shutil.move(source_path, destination_path)
    logging.info(f"Successfully moved {source_path} to {destination_path}.")

def process_single_ligand(args_tuple):
    """Worker function to process one ligand. Designed for multiprocessing."""
    ligand_name, pdb, work_path, docking_path = args_tuple

    # Each process should log to its own file to avoid garbled output
    log_file_path = os.path.join(docking_path, ligand_name, f"{ligand_name}_processing.log")

    try:
        logging.info(f"--- Starting processing for ligand: {ligand_name} ---")

        # Setup directories and files for this specific ligand
        ligand_docking_dir = os.path.join(docking_path, ligand_name)
        if os.path.exists(ligand_docking_dir):
            shutil.rmtree(ligand_docking_dir)
        os.makedirs(ligand_docking_dir)

        # We need the original prepared SDF to extract this one molecule
        all_ligands_sdf_path = os.path.join(work_path, "uff", f"{pdb}_uff.sdf")
        single_sdf_path = os.path.join(ligand_docking_dir, f"{ligand_name}.sdf")

        found = False
        with Chem.SDMolSupplier(all_ligands_sdf_path) as supplier:
            for mol in supplier:
                if mol and mol.HasProp("_Name") and mol.GetProp("_Name") == ligand_name:
                    with Chem.SDWriter(single_sdf_path) as writer:
                        writer.write(mol)
                    found = True
                    break
        if not found:
            raise PipelineStepFailed(f"Could not find molecule for '{ligand_name}' in SDF.")

        # Copy necessary base files/folders
        for required_dir in ["ligand", "pocket"]:
            shutil.copytree(os.path.join(work_path, required_dir), os.path.join(ligand_docking_dir, required_dir))

        shutil.copytree(os.path.join(work_path, "uff"), os.path.join(ligand_docking_dir, "uff"))

        # 1. Create query CSV
        single_uff_infer_csv = os.path.join(ligand_docking_dir, f"{ligand_name}_infer.csv")
        run_command(["python", "tools/inference/inter_sdf2csv.py", single_sdf_path, "1"], log_file_path)
        generated_csv_path = single_sdf_path.replace('.sdf', '_infer.csv')
        if os.path.exists(generated_csv_path):
            os.replace(generated_csv_path, single_uff_infer_csv)

        # 2. Predict energy
        run_command(["python", "inference.py", "-test_csv", single_uff_infer_csv, "-work_path", ligand_docking_dir, "-ensemble", "../checkpoints/v0.2_energy_model", "-gpus", "0", "-batch_size", "1", "-energy_output_folder", ".", "-uff_as_ligand"], log_file_path)

        # 3. Reconstruct (Docking)
        run_command(["python", "docking/reconstruct_ligands.py", "-y", "--cwd", ligand_docking_dir, "-i", pdb, "--uff_folder", "uff", "find"], log_file_path)

        # 4. Scoring
        reconstructed_sdf = os.path.join(ligand_docking_dir, "ligand_reconstructing", f"{pdb}_docked.sdf")
        if not os.path.exists(reconstructed_sdf):
             raise PipelineStepFailed(f"Reconstructed SDF not found for {ligand_name}")

        run_command(["python", "tools/inference/inter_sdf2csv.py", reconstructed_sdf, "0"], log_file_path)
        docked_infer_csv = reconstructed_sdf.replace('.sdf', '_infer.csv')

        # The scoring script saves to a hardcoded `result` dir. We handle this in the main thread.
        run_command(["python", "inference.py", "-test_csv", docked_infer_csv, "-work_path", work_path, "-ligand_folder", os.path.join(os.path.relpath(ligand_docking_dir, work_path), "ligand_reconstructing"), "-ensemble", "checkpoints/v0.2_affinity_model/model*", "-vs", "-gpus", "0", "-batch_size", "20", "--pose_sel", "True"], log_file_path)

        # Move the final result to a unique location
        final_result_csv = os.path.join("result", f"{pdb}_docked_infer_ensemble.csv")
        unique_result_csv = os.path.join(ligand_docking_dir, f"{ligand_name}_result.csv")
        if os.path.exists(final_result_csv):
            os.replace(final_result_csv, unique_result_csv)

        return (ligand_name, "success", None)

    except Exception as e:
        logging.error(f"Caught exception for ligand {ligand_name}: {e}")
        return (ligand_name, "failed", str(e))

def run_docking_parallel(work_path, docking_path, pdb, omp_num_threads, ligands_to_process, num_workers):
    """Runs the docking process in parallel for a list of ligands."""
    logging.info(f"Starting parallel docking with {num_workers} workers.")

    # Prepare arguments for each worker
    tasks = [(ligand, pdb, work_path, docking_path) for ligand in ligands_to_process]

    docking_success = []
    docking_failed = {}

    with Pool(processes=num_workers) as pool:
        results = pool.map(process_single_ligand, tasks)

    for res in results:
        ligand_name, status, error_msg = res
        if status == "success":
            docking_success.append(ligand_name)
        else:
            docking_failed[ligand_name] = error_msg

    # Merge results
    logging.info("--- Merging individual docking results ---")
    all_dfs = []
    for ligand_name in docking_success:
        unique_result_csv = os.path.join(docking_path, ligand_name, f"{ligand_name}_result.csv")
        if os.path.exists(unique_result_csv):
            df = pd.read_csv(unique_result_csv)
            df['ligand_name'] = ligand_name # Ensure we can track origin
            all_dfs.append(df)

    if all_dfs:
        merged_df = pd.concat(all_dfs, ignore_index=True)
        merged_csv_path = os.path.join("result", f"{pdb}_docked_infer_ensemble.csv")
        merged_df.to_csv(merged_csv_path, index=False)
        logging.info(f"Merged ensemble CSV saved to {merged_csv_path}")

    final_docked_sdf_path = os.path.join(work_path, "infer", f"{pdb}_docked.sdf")
    with Chem.SDWriter(final_docked_sdf_path) as writer:
        for ligand_name in docking_success:
            reconstructed_sdf = os.path.join(docking_path, ligand_name, "ligand_reconstructing", f"{pdb}_docked.sdf")
            if os.path.exists(reconstructed_sdf):
                for mol in Chem.SDMolSupplier(reconstructed_sdf):
                    if mol: writer.write(mol)
    logging.info(f"Merged docked SDF saved to {final_docked_sdf_path}")

    return docking_success, docking_failed

def analyze_results(work_path, pdb, df_lig):
    """Analyzes the docking results."""
    logging.info("Analyzing results...")
    results_csv = os.path.join("result", f"{pdb}_docked_infer_ensemble.csv")
    if not os.path.exists(results_csv):
        logging.warning(f"Results CSV not found: {results_csv}. Skipping analysis.")
        return
    log_path = os.path.join("result", "analysis_commands.log")
    run_command(["python", "analyze.py", "--results-csv", results_csv, "--original-csv", os.path.join(work_path, "raw", f"{df_lig}.csv"), "--experimental-col", "pValue", "--output-folder", "result"], log_path)

def main():
    parser = argparse.ArgumentParser(description="Master pipeline for Interformer docking with resume and parallel capability.")
    parser.add_argument("-p", "--pdb", required=True, help="PDB code.")
    parser.add_argument("-w", "--work_path", required=True, help="Working directory.")
    parser.add_argument("-d", "--docking_path", required=True, help="Directory for docking results.")
    parser.add_argument("-l", "--df_lig", required=True, help="Name of the CSV file with ligands.")
    parser.add_argument("-o", "--omp_num_threads", default="1,1", help="Number of threads for OpenMP (less relevant in parallel mode).")
    parser.add_argument("--parallel_workers", type=int, default=0, help="Number of parallel workers. 0 means use all available CPUs minus one.")
    args = parser.parse_args()

    num_workers = args.parallel_workers if args.parallel_workers > 0 else max(1, cpu_count() - 1)
    logging.info(f"Starting master pipeline. Using {num_workers} parallel workers.")

    report_csv_path = os.path.join("result", f"pipeline_summary_report_{args.pdb}.csv")

    # Initial setup is always sequential
    try:
        prepare_directories(args.work_path)
        copy_input_files(args.work_path, args.pdb, args.df_lig)
        prepare_protein_and_ligand(args.work_path, args.pdb)
    except (PipelineStepFailed, IOError, OSError) as e:
        logging.critical(f"Setup phase failed: {e}. Aborting.")
        sys.exit(1)

    successful_prep, failed_prep = prepare_ligands_for_interformer(os.path.join(args.work_path, "raw", f"{args.df_lig}.csv"), os.path.join(args.work_path, "uff", f"{args.pdb}_uff.sdf"))

    ligands_to_process = successful_prep

    if os.path.exists(report_csv_path):
        logging.info(f"Found existing report file: {report_csv_path}. Attempting to resume.")
        report_df = pd.read_csv(report_csv_path)
        pipeline_report = report_df.to_dict('records')
        completed_ligands = set(report_df[report_df['docking'].isin(['success', 'failed'])]['ligand_name'])
        ligands_to_process = [lig for lig in successful_prep if lig not in completed_ligands]
        logging.info(f"Skipping {len(completed_ligands)} already processed ligands. Processing {len(ligands_to_process)} new ligands.")
    else:
        pipeline_report = []
        for ligand in successful_prep:
            pipeline_report.append({"ligand_name": ligand, "preparation": "success", "docking": "pending", "analysis": "pending", "error": ""})
        for ligand in failed_prep:
            pipeline_report.append({"ligand_name": ligand, "preparation": "failed", "docking": "skipped", "analysis": "skipped", "error": "Failed during ligand preparation"})

    if ligands_to_process:
        docking_success, docking_failed = run_docking_parallel(args.work_path, args.docking_path, args.pdb, args.omp_num_threads, ligands_to_process, num_workers)

        report_map = {item['ligand_name']: item for item in pipeline_report}
        for ligand_name in docking_success:
            if ligand_name in report_map: report_map[ligand_name]['docking'] = 'success'
        for ligand_name, error_msg in docking_failed.items():
            if ligand_name in report_map:
                report_map[ligand_name].update({'docking': 'failed', 'analysis': 'skipped', 'error': error_msg})

        try:
            if docking_success:
                analyze_results(args.work_path, args.pdb, args.df_lig)
                for ligand_name in docking_success:
                    if ligand_name in report_map: report_map[ligand_name]['analysis'] = 'success'
        except PipelineStepFailed as e:
            logging.error(f"Analysis phase failed: {e}")
            for ligand_name in docking_success:
                if ligand_name in report_map: report_map[ligand_name].update({'analysis': 'failed', 'error': str(e)})

    final_report_df = pd.DataFrame(pipeline_report)
    final_report_df.to_csv(report_csv_path, index=False)
    logging.info(f"Pipeline summary report saved to {report_csv_path}")
    logging.info("Master pipeline finished.")

if __name__ == "__main__":
    main()
