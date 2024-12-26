import os
import subprocess
import argparse
import pandas as pd
import re
from Bio import SeqIO

def run_command(command):
    """Run a shell command and return the output, error, and exit code."""
    try:
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return result.stdout, result.stderr, result.returncode
    except Exception as e:
        print(f"Exception occurred while executing command: {command}\n{e}")
        return None, None, 1

def index_reference(reference_path):
    """Index the reference file using BWA if not already indexed."""
    index_files = [reference_path + ext for ext in ['.bwt', '.pac', '.ann', '.amb', '.sa']]
    
    if all(os.path.exists(idx) for idx in index_files):
        print(f"Reference {reference_path} is already indexed.")
        return True
        
    command = f"bwa index {reference_path}"
    print(f"Indexing {reference_path}...")
    _, stderr, code = run_command(command)
    
    if code == 0:
        print(f"Indexing complete for {reference_path}")
        return True
    else:
        print(f"Indexing failed for {reference_path}: {stderr}")
        return False

def calculate_correctness(aligned_sequence, output_sequence):
    """Calculate correctness as the percentage of matching characters (ignoring gaps and case)."""
    aligned_sequence = re.sub(r"[-\s]", "", aligned_sequence).upper()
    output_sequence = re.sub(r"[-\s]", "", output_sequence).upper()
    
    if not aligned_sequence or not output_sequence:
        return 0.0
    
    matches = sum(1 for a, b in zip(aligned_sequence, output_sequence) if a == b)
    correctness = (matches / len(aligned_sequence)) * 100 if aligned_sequence else 0
    return correctness

def parse_alignment_output(stdout, stderr, query_name):
    """Parse alignment output from BWA logs."""
    execution_time = None

    try:
        for line in stderr.split("\n"):
            if "[main] Real time" in line:
                execution_time = float(re.search(r"Real time:\s+([\d.]+)", line).group(1))

    except Exception as e:
        print(f"Error parsing alignment output: {e}")

    return {
        "Sequence_Name": query_name[4:-6],
        "Execution_Time (sec)": execution_time,
    }

def align_method(query_path, country, fasta_file, method, output_dir, aligned_folder, reference_path):
    """Run a BWA alignment method (MEM or SW) and process results."""
    # Get sequence length first
    seq_length = None
    with open(query_path, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            seq_length = len(record.seq)
            break

    if method == "MEM":
        command = f"bwa mem {reference_path} {query_path}"
    elif method == "SW":
        command = f"bwa bwasw {reference_path} {query_path}"
    else:
        raise ValueError(f"Invalid method: {method}")

    print(f"Running BWA {method} for {query_path}...")
    stdout, stderr, returncode = run_command(command)

    if returncode != 0:
        print(f"Error running BWA {method} for {query_path}: {stderr}")
        return None

    # Extract the aligned sequence from the SAM output
    output_sequence = None
    for line in stdout.split("\n"):
        if not line.startswith("@") and line.strip():
            fields = line.split("\t")
            output_sequence = fields[9]  # The aligned sequence is in the 10th column
            break

    # Read the corresponding aligned sequence
    aligned_path = os.path.join(aligned_folder, fasta_file)
    aligned_sequence = None
    if os.path.exists(aligned_path):
        for record in SeqIO.parse(aligned_path, "fasta"):
            aligned_sequence = str(record.seq)
            break

    # Calculate correctness metric
    correctness = calculate_correctness(aligned_sequence, output_sequence) if aligned_sequence else None

    result = parse_alignment_output(stdout, stderr, fasta_file)
    result["Country"] = country.capitalize()
    result["Correctness (%)"] = correctness
    result["Sequence_Length"] = seq_length  # Add sequence length to results
    return result

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
        
        # Run both MEM and SW methods
        for method in ["MEM", "SW"]:
            result = align_method(unaligned_file, country_name, fasta_file, 
                                method, output_dir, aligned_folder, reference_path)
            if result:
                results_dict[method].append(result)

def main():
    parser = argparse.ArgumentParser(description="Run BWA alignment pipeline using specified methods.")
    parser.add_argument('-i', '--input_dir', type=str, required=True, 
                       help='Path to the root directory containing country folders.')
    parser.add_argument('-o', '--output_dir', type=str, default="results", 
                       help='Path to the directory where CSV files will be saved.')
    parser.add_argument('-r', '--reference', type=str, default='reference/reference.fasta', 
                       help='Path to the reference FASTA file.')
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir
    reference_path = args.reference

    os.makedirs(output_dir, exist_ok=True)

    # Initialize results dictionary for each algorithm
    results_dict = {
        "MEM": [],
        "SW": []
    }

    for country in os.listdir(input_dir):
        country_path = os.path.join(input_dir, country)
        if os.path.isdir(country_path):
            print(f"Processing country: {country}")
            process_country_folder(country_path, results_dict, output_dir, reference_path)

    column_order = [
        "Sequence_Name",
        "Country",
        "Sequence_Length",
        "Execution_Time (sec)",
        "Correctness (%)"
    ]

    # Save results to separate CSV files with ordered columns
    for method, results in results_dict.items():
        csv_path = os.path.join(output_dir, f"{method}_results.csv")
        df = pd.DataFrame(results)
        # Reorder columns to match backtrack CSV
        df = df[column_order]
        df.to_csv(csv_path, index=False)
        print(f"\n{method} results have been saved to {csv_path}")

if __name__ == "__main__":
    main()