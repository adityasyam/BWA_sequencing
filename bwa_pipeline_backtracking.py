import os
import subprocess
import argparse
import pandas as pd
import re
from Bio import SeqIO
import time

def run_command_with_timing(command):
    """Run a shell command and return the output, error, exit code, and execution time."""
    start_time = time.time()
    try:
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        end_time = time.time()
        execution_time = end_time - start_time
        return result.stdout, result.stderr, result.returncode, execution_time
    except subprocess.CalledProcessError as e:
        print(f"Command failed with return code {e.returncode}: {command}\n{e.stderr}")
        return None, e.stderr, e.returncode, 0
    except Exception as e:
        print(f"Exception occurred while executing command: {command}\n{e}")
        return None, str(e), 1, 0

def index_reference(reference_path='reference/reference.fasta'):
    """Index the reference file using BWA if not already indexed."""
    index_files = [reference_path + ext for ext in ['.bwt', '.pac', '.ann', '.amb', '.sa']]

    if all(os.path.exists(idx) for idx in index_files):
        print(f"Reference {reference_path} is already indexed.")
        return True

    command = f"bwa index {reference_path}"
    print(f"Indexing {reference_path}...")
    _, stderr, code, _ = run_command_with_timing(command)

    if code == 0:
        print(f"Indexing complete for {reference_path}")
        return True
    else:
        print(f"Indexing failed for {reference_path}: {stderr}")
        return False

def calculate_correctness(aligned_sequence, output_sequence):
    """Calculate correctness as the percentage of matching characters (ignoring gaps and case)."""
    if not aligned_sequence or not output_sequence:
        return 0.0
    
    aligned_sequence = re.sub(r"[-\s]", "", aligned_sequence).upper()
    output_sequence = re.sub(r"[-\s]", "", output_sequence).upper()
    
    matches = sum(1 for a, b in zip(aligned_sequence, output_sequence) if a == b)
    correctness = (matches / len(aligned_sequence)) * 100 if aligned_sequence else 0
    return correctness

def align_backtrack(query_path, country, fasta_file, output_dir, aligned_folder, reference_path):
    """Run BWA Backtrack alignment with comprehensive timing."""
    try:
        sai_file = os.path.join(output_dir, f"{fasta_file}.sai")
        
        # Time the aln step
        aln_command = f"bwa aln -t 1 {reference_path} {query_path} > {sai_file}"
        print(f"Running BWA aln for {query_path}...")
        _, stderr_aln, returncode_aln, aln_time = run_command_with_timing(aln_command)

        if returncode_aln != 0:
            print(f"Error running BWA aln for {query_path}: {stderr_aln}")
            return None

        # Time the samse step
        samse_command = f"bwa samse {reference_path} {sai_file} {query_path}"
        print(f"Running BWA samse for {query_path}...")
        stdout_samse, stderr_samse, returncode_samse, samse_time = run_command_with_timing(samse_command)

        if returncode_samse != 0:
            print(f"Error running BWA samse for {query_path}: {stderr_samse}")
            return None

        # Extract sequence from SAM output
        output_sequence = None
        for line in stdout_samse.split("\n"):
            if not line.startswith("@") and line.strip():
                fields = line.split("\t")
                if len(fields) >= 10:
                    output_sequence = fields[9]
                    break

        # Read aligned sequence
        aligned_path = os.path.join(aligned_folder, fasta_file)
        aligned_sequence = None
        if os.path.exists(aligned_path):
            for record in SeqIO.parse(aligned_path, "fasta"):
                aligned_sequence = str(record.seq)
                break

        # Calculate correctness
        correctness = calculate_correctness(aligned_sequence, output_sequence) if aligned_sequence and output_sequence else None

        # Calculate and round total time
        total_time = round(aln_time + samse_time, 3)

        # Get sequence length for verification
        with open(query_path, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                seq_length = len(record.seq)
                break

        return {
            "Sequence_Name": fasta_file[4:-6],
            "Country": country.capitalize(),
            "Sequence_Length": seq_length,
            "Execution_Time (sec)": total_time,
            "Correctness (%)": correctness
        }
    
    except Exception as e:
        print(f"Error in align_backtrack for {query_path}: {str(e)}")
        return None

def process_country_folder(country_path, results_dict, output_dir, reference_path):
    """Process a country's folder and align all paired FASTA files."""
    unaligned_folder = os.path.join(country_path, "unaligned_seq")
    aligned_folder = os.path.join(country_path, "aligned_seq")
    country_name = os.path.basename(country_path)

    if not os.path.isdir(unaligned_folder):
        print(f"Skipping {country_path}: No unaligned_seq folder found.")
        return

    if not os.path.isdir(aligned_folder):
        print(f"Skipping {country_path}: No aligned_seq folder found.")
        return

    if not index_reference(reference_path):
        print(f"Could not index reference {reference_path}. Skipping country.")
        return

    fasta_files = [f for f in os.listdir(unaligned_folder) if f.endswith('.fasta')]

    for fasta_file in fasta_files:
        unaligned_file = os.path.join(unaligned_folder, fasta_file)
        backtrack_result = align_backtrack(unaligned_file, country_name, fasta_file, output_dir, aligned_folder, reference_path)
        if backtrack_result:
            results_dict["Backtrack"].append(backtrack_result)

def main():
    parser = argparse.ArgumentParser(description="Run BWA alignment pipeline using a specified reference.")
    parser.add_argument('-i', '--input_dir', type=str, required=True, help='Path to the root directory containing country folders.')
    parser.add_argument('-o', '--output_dir', type=str, default="results", help='Path to the directory where intermediate files and CSV will be saved.')
    parser.add_argument('-r', '--reference', type=str, default='reference/reference.fasta', help='Path to the reference FASTA file.')
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir
    reference_path = args.reference

    os.makedirs(output_dir, exist_ok=True)

    results_dict = {"Backtrack": []}

    for country in os.listdir(input_dir):
        country_path = os.path.join(input_dir, country)
        if os.path.isdir(country_path):
            print(f"Processing country: {country}")
            process_country_folder(country_path, results_dict, output_dir, reference_path)

    # Create DataFrame with all timing information
    csv_path = os.path.join(output_dir, "Backtrack_results.csv")
    df = pd.DataFrame(results_dict["Backtrack"])
    df.to_csv(csv_path, index=False)
    print(f"\nBacktrack results have been saved to {csv_path}")

if __name__ == "__main__":
    main()