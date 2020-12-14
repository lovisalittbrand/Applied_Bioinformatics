PIPELINE

Within this folder are all necessary files to run the pipeline located. The pipeline aims of performing the following set of analysis on arCOGS:
- Alignment (MAFFT)
- Trimming (TrimAl)
- Concatenation + mapping (upp_concatArCOGalignments.py)
- Phylogenetic reconstruction (Iqtree)
- Frequency + Entropy visualization (align_vis.py) and phylogenetic signal visualization (signal_vis2trees.py)

The pipeline is written in bash and can be found in bash_pipeline.sh, where it is clearly explained how the pipeline should be run. The pipeline can be run in two ways which is exemplified in run_bash_p1.sh (part 1) and run_bash_p2.sh (part 2). Part 1 should be run in case you have unaligned FASTA-files of arCOGS. Part 2 should be run if you have a concatenation of arCOGs and a tree-file containing a topology constraint that will be used in the run for Iqtree. The result from both parts will be handled in run_bash_p2.sh where a file containing site likelihoods in the correct format is created and thereafter visualized. 
