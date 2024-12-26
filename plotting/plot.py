import pandas as pd
import matplotlib.pyplot as plt

def plot_execution_time(csv_path_backtracking, csv_path_mem, csv_path_sw, save_path, plot_title):

    # Load the CSV files
    backtrack_results = pd.read_csv(csv_path_backtracking)
    mem_results = pd.read_csv(csv_path_mem)
    sw_results = pd.read_csv(csv_path_sw)

    # Merge the dataframes and annotate with the variant
    backtrack_results['Variant'] = 'Backtracking'
    mem_results['Variant'] = 'MEM'
    sw_results['Variant'] = 'SW'

    # Concatenate all results into a single dataframe for easier plotting
    combined_results = pd.concat([backtrack_results, mem_results, sw_results])

    # Create a line plot
    plt.figure(figsize=(12, 8))
    for variant, group in combined_results.groupby('Variant'):
        # Group by sequence name and calculate mean of numeric columns
        grouped_by_sequence = group.groupby('Sequence_Name').mean(numeric_only=True)
        plt.plot(
            grouped_by_sequence.index,
            grouped_by_sequence['Execution_Time (sec)'],
            label=f"{variant}",
            marker='o',
            alpha=0.5 if variant == 'SW' else 1.0  # Fixed opacity for SW variant
        )

    # Configure axis labels and title
    plt.xticks(rotation=45, ha='right', fontsize=9)
    plt.xlabel('Sequences', fontsize=12)
    plt.ylabel('Execution Time (s)', fontsize=12)
    plt.title(plot_title, fontsize=14)
    plt.legend(title="Variants", fontsize=10)

    # Save and show the plot
    plt.tight_layout()
    plt.savefig(save_path)


# Example usage
# REMOVE COMMENTS ON THE DESIRED BASE PAIR LENGTH TO RUN THE CORRESPONDING PLOTTING

# plot_execution_time(
#     '../results/alignment_result_50.csv/Backtrack_results.csv',
#     '../results/alignment_result_50_ms.csv/MEM_results.csv',
#     '../results/alignment_result_50_ms.csv/SW_results.csv',
#     'graphs/50bp_results.png',
#     'Execution Time for BWA Variants for 50bp'
# )

# plot_execution_time(
#     '../results/alignment_result_200.csv/Backtrack_results.csv',
#     '../results/alignment_result_200_ms.csv/MEM_results.csv',
#     '../results/alignment_result_200_ms.csv/SW_results.csv',
#     'graphs/200bp_results.png',
#     'Execution Time for BWA Variants for 200bp'
# )

# plot_execution_time(
#     '../results/alignment_result_1k.csv/Backtrack_results.csv',
#     '../results/alignment_result_1k_ms.csv/MEM_results.csv',
#     '../results/alignment_result_1k_ms.csv/SW_results.csv',
#     'graphs/1000bp_results.png',
#     'Execution Time for BWA Variants for 1000bp'
# )

# plot_execution_time(
#     '../results/alignment_result_2k.csv/Backtrack_results.csv',
#     '../results/alignment_result_2k_ms.csv/MEM_results.csv',
#     '../results/alignment_result_2k_ms.csv/SW_results.csv',
#     'graphs/2000bp_results.png',
#     'Execution Time for BWA Variants for 2000bp'
# )

# plot_execution_time(
#     '../results/alignment_result_5k.csv/Backtrack_results.csv',
#     '../results/alignment_result_5k_ms.csv/MEM_results.csv',
#     '../results/alignment_result_5k_ms.csv/SW_results.csv',
#     'graphs/5000bp_results.png',
#     'Execution Time for BWA Variants for 5000bp'
# )

# plot_execution_time(
#     '../results/alignment_result_9k.csv/Backtrack_results.csv',
#     '../results/alignment_result_9k_ms.csv/MEM_results.csv',
#     '../results/alignment_result_9k_ms.csv/SW_results.csv',
#     'graphs/9000bp_results.png',
#     'Execution Time for BWA Variants for 9000bp'
# )

# plot_execution_time(
#     '../results/alignment_result_9749.csv/Backtrack_results.csv',
#     '../results/alignment_result_9749_ms.csv/MEM_results.csv',
#     '../results/alignment_result_9749_ms.csv/SW_results.csv',
#     'graphs/9749bp_results.png',
#     'Execution Time for BWA Variants for 9749bp'
# )
