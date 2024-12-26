This is a code repository to run and benchmark different versions of BWA for DNA sequencing, particularly on HIV-1 nucleotide sequences.

Prerequisites:
Create a virtual environment and install the requirements listed in the `requirements.txt` file using
`pip install -r requirements.txt`

Then, in the same virtual environment, run
`brew install bwa`

The current genome sequences used are as follows:
The reference sequence is saved in `reference/reference.fasta`
The current reference being used is the complete HXB2 genome sequence from LANL. 
To download the reference sequence, navigate to `https://www.hiv.lanl.gov/components/sequence/HIV/asearch/query_one.comp?se_id=K03455`, check "Entire Sequence" and click "Download"


To obtain the query sequences, first navigate to `https://www.hiv.lanl.gov/content/index` and click on "Sequence Database". Then, follow these steps:
1) Use the "Search Interface" option 
2) This will take you to the "Sequence Search Interface" page
3) Enter the desired "Sequence length" (starting from 50 and going up to 9749)
4) Ensure "Virus" tab says "HIV-1"
5) Click "Search"
6) Select all the sequences you want to download, then click on "Download sequences"
7) Ensure the "Sequence type" says "Nucleotides"
8) Choose the "Fasta" option under the "Format" button
9) Ensure "Pad aligned sequences" is unchecked 
10) Ensure "Gap handling" is set to "none"
11) Make sure the "Align" option is unchecked when you download unaligned (query) sequences
12) Make sure the "Align" option is checked when you download aligned (ground truth) sequences
13) Click "Download" 
14) You will be taken to a fresh page with a fasta file linked-click on the fasta link and the file will be downloaded to your device
15) The downloaded fasta file contains all the sequence you selected, in order of their display on the LANL page

There are different data folders for different sequence lengths:
1) `data/data_50` for sequence length 50
2) `data/data_200` for sequence length 200
3) `data/data_1k` for sequence length 1000
4) `data/data_2k` for sequence length 2000
5) `data/data_5k` for sequence length 5000
6) `data/data_9k` for sequence length 9000
7) `data/data_9749` for sequence length 9749 (the longest sequence found in the database)

For each geographic region (country), we constructed two datasets: a set of unaligned query sequences and their corresponding unpadded aligned sequences which served as our ground truth dataset to verify correctness. We saved the unaligned fasta file in a directory adjacent to the directory containing the aligned fasta files, allowing for parallel verification later in our data pipeline. Essentially, we had identically named fasta files for each query sequence in the `aligned_seq` and `unaligned_seq` directories for each country. The Assumed Directory Structure Sample at the top clearly shows this.


**If you wish to construct datasets with different data, you will need to replicate a similar directory structure to streamline the sequencing pipeline.

Execution:
To run the BWA algorithms, choose one of the data folders corresponding to a sequence length and named the output CSV file as desired.

For backtracking run the following command:
`python bwa_pipeline_backtracking.py -i data/path_to_specific_folder -o results/path_to_specific_alignment_result.csv`

For MEM/SW run this command:
`python bwa_pipeline_mem_sw.py -i data/path_to_specific_folder -o results/path_to_specific_alignment_result.csv`

There are different output CSV folders for different sequence lengths:
1) `results/alignment_result_50.csv` for backtracking alignment for sequence length 50
2) `results/alignment_result_50_ms.csv` for MEM/SW alignment for sequence length 50
3) `results/alignment_result_200.csv` for backtracking alignment for sequence length 200
4) `results/alignment_result_200_ms.csv` for MEM/SW alignment for sequence length 200
5) `results/alignment_result_1k.csv` for backtracking alignment for sequence length 1000
6) `results/alignment_result_1k_ms.csv` for MEM/SW alignment for sequence length 1000
7) `results/alignment_result_2k.csv` for backtracking alignment for sequence length 2000
8) `results/alignment_result_2k_ms.csv` for MEM/SW alignment for sequence length 2000
9) `results/alignment_result_5k.csv` for backtracking alignment for sequence length 5000
10) `results/alignment_result_5k_ms.csv` for MEM/SW alignment for sequence length 5000
11) `results/alignment_result_9k.csv` for backtracking alignment for sequence length 9000
12) `results/alignment_result_9k_ms.csv` for MEM/SW alignment for sequence length 9000
13) `results/alignment_result_9749.csv` for backtracking alignment for sequence length 9749
14) `results/alignment_result_9749_ms.csv` for MEM/SW alignment for sequence length 9749

Plotting:
To recreate the plots in `plotting/graphs` after any changes in the output CSV files,
first uncomment the section of code in `plotting/plot.py` that corresponds to the desired sequence length and then run
`cd plotting`
`python plot.py`

The appropriate image in the `plotting/graphs` directory will be automatically overwritten.

Statistical Tests:
Finally, our repository also contains a script to compute the p-values and significance of the performance metrics. This is written in `tests/test.py`.
To run the statistical test script, first uncomment the section of code in the script that corresponds to the desired sequence length and then run
`cd tests`
`python test.py`

The resultant CSV files are saved in the adjacent `tests/test_outputs` directory and will be automatically overwritten if rerun. The current output files are as follows:
1) `test_outputs/50bp_tests.csv` for statistical analysis for 50bp
2) `test_outputs/200bp_tests.csv` for statistical analysis for 200bp
3) `test_outputs/1000bp_tests.csv` for statistical analysis for 1000bp
4) `test_outputs/2000bp_tests.csv` for statistical analysis for 2000bp
5) `test_outputs/5000bp_tests.csv` for statistical analysis for 5000bp
6) `test_outputs/9000bp_tests.csv` for statistical analysis for 9000bp
7) `test_outputs/9749bp_tests.csv` for statistical analysis for 9749bp

**References**:
1) Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics, 25, 1754-1760.

2) Li H. and Durbin R. (2010) Fast and accurate long-read alignment with Burrows-Wheeler transform. Bioinformatics, 26, 589-595.

3) Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v2 [q-bio.GN].
