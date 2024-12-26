import pandas as pd
import numpy as np
from scipy import stats
import argparse

def load_and_prepare_data(backtrack_path, mem_path, sw_path):
    """Load data from CSV files and prepare for analysis."""
    backtrack_df = pd.read_csv(backtrack_path)
    mem_df = pd.read_csv(mem_path)
    sw_df = pd.read_csv(sw_path)
    
    # Add algorithm identifier
    backtrack_df['Algorithm'] = 'Backtrack'
    mem_df['Algorithm'] = 'MEM'
    sw_df['Algorithm'] = 'SW'
    
    return backtrack_df, mem_df, sw_df

def perform_statistical_tests(backtrack_df, mem_df, sw_df):
    """Perform statistical tests between algorithms."""
    results = []
    
    # Metrics to test - modified names
    metrics = ['Execution_Time (sec)', 'Correctness (%)']
    display_names = {'Execution_Time (sec)': 'Time', 'Correctness (%)': 'Correctness'}
    
    # Algorithm pairs to compare
    algorithm_pairs = [
        ('Backtrack', 'MEM'),
        ('Backtrack', 'SW'),
        ('MEM', 'SW')
    ]
    
    for metric in metrics:
        for alg1, alg2 in algorithm_pairs:
            # Get data for each algorithm
            data1 = backtrack_df[metric] if alg1 == 'Backtrack' else \
                   mem_df[metric] if alg1 == 'MEM' else sw_df[metric]
            data2 = backtrack_df[metric] if alg2 == 'Backtrack' else \
                   mem_df[metric] if alg2 == 'MEM' else sw_df[metric]
            
            # Perform statistical tests
            # Mann-Whitney U test for non-parametric comparison
            statistic, p_value = stats.mannwhitneyu(data1, data2, alternative='two-sided')
            
            # Calculate effect size (Cohen's d)
            cohens_d = (np.mean(data1) - np.mean(data2)) / \
                      np.sqrt((np.std(data1) ** 2 + np.std(data2) ** 2) / 2)
            
            results.append({
                'Metric': display_names[metric],  # Use simplified metric name
                'Algorithm_1': alg1,
                'Algorithm_2': alg2,
                'Mean_1': round(np.mean(data1), 4),
                'Mean_2': round(np.mean(data2), 4),
                'P_Value': round(p_value, 4),
                'Effect_Size': round(cohens_d, 4),
                'Significant': p_value < 0.05
            })
    
    return pd.DataFrame(results)

def add_interpretation(results_df):
    """Add interpretation column to results."""
    def interpret_result(row):
        if not row['Significant']:
            return "N/A"
        
        effect_size = abs(row['Effect_Size'])
        if effect_size < 0.2:
            magnitude = "negligible"
        elif effect_size < 0.5:
            magnitude = "small"
        elif effect_size < 0.8:
            magnitude = "medium"
        else:
            magnitude = "large"
            
        better_alg = row['Algorithm_1'] if row['Mean_1'] > row['Mean_2'] else row['Algorithm_2']
        worse_alg = row['Algorithm_2'] if row['Mean_1'] > row['Mean_2'] else row['Algorithm_1']
        
        if row['Metric'] == 'Time':  # Updated metric name in condition
            better_alg, worse_alg = worse_alg, better_alg
            
        return f"{better_alg}>{worse_alg}, {magnitude}"
    
    results_df['Interpretation'] = results_df.apply(interpret_result, axis=1)
    return results_df


def main():

    # backtrack_path='../results/alignment_result_50.csv/Backtrack_results.csv'
    # mem_path='../results/alignment_result_50_ms.csv/MEM_results.csv'
    # sw_path='../results/alignment_result_50_ms.csv/SW_results.csv'
    # output_path='test_outputs/50bp_tests.csv'

    # backtrack_path='../results/alignment_result_200.csv/Backtrack_results.csv'
    # mem_path='../results/alignment_result_200_ms.csv/MEM_results.csv'
    # sw_path='../results/alignment_result_200_ms.csv/SW_results.csv'
    # output_path='test_outputs/200bp_tests.csv'

    # backtrack_path='../results/alignment_result_1k.csv/Backtrack_results.csv'
    # mem_path='../results/alignment_result_1k_ms.csv/MEM_results.csv'
    # sw_path='../results/alignment_result_1k_ms.csv/SW_results.csv'
    # output_path='test_outputs/1000bp_tests.csv'

    # backtrack_path='../results/alignment_result_2k.csv/Backtrack_results.csv'
    # mem_path='../results/alignment_result_2k_ms.csv/MEM_results.csv'
    # sw_path='../results/alignment_result_2k_ms.csv/SW_results.csv'
    # output_path='test_outputs/2000bp_tests.csv'

    # backtrack_path='../results/alignment_result_5k.csv/Backtrack_results.csv'
    # mem_path='../results/alignment_result_5k_ms.csv/MEM_results.csv'
    # sw_path='../results/alignment_result_5k_ms.csv/SW_results.csv'
    # output_path='test_outputs/5000bp_tests.csv'

    # backtrack_path='../results/alignment_result_9k.csv/Backtrack_results.csv'
    # mem_path='../results/alignment_result_9k_ms.csv/MEM_results.csv'
    # sw_path='../results/alignment_result_9k_ms.csv/SW_results.csv'
    # output_path='test_outputs/9000bp_tests.csv'

    # backtrack_path='../results/alignment_result_9749.csv/Backtrack_results.csv'
    # mem_path='../results/alignment_result_9749_ms.csv/MEM_results.csv'
    # sw_path='../results/alignment_result_9749_ms.csv/SW_results.csv'
    # output_path='test_outputs/9749bp_tests.csv'

    backtrack_df, mem_df, sw_df = load_and_prepare_data(
        backtrack_path, mem_path, sw_path)
    
    # Perform statistical analysis
    results_df = perform_statistical_tests(backtrack_df, mem_df, sw_df)
    
    # Add interpretations
    results_df = add_interpretation(results_df)
    
    # Save results
    results_df.to_csv(output_path, index=False)
    print(f"Analysis results saved to {output_path}")

if __name__ == "__main__":
    main()