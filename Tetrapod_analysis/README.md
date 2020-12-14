# TETRAPOD ANALYSIS

This folder contains the raw data from the tetrapod analysis (0_raw_data) and the results from the analysis in the other two folders. The folders contain the following:

- 0_raw_data: raw data of tetrapods in format of unaligned FASTA-files of a set of genes
- 1_untrimmed: Result after aligning the raw data sequences with MAFFT without performing any trimming. The files are thereafter concatenated and can be found in original_concatenation-folder. Thereafter, do the other folders contain the result from the analysis after modifying the original concatenation in several ways.
- 2_trimmed: Result after aligning the raw data sequences as well as performing trimming with TrimAl (gappyout). The files are therafter concatenated. However, no further analysis has yet been done on the concatenation.
