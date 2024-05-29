# Type-2-diabetes-mellitus-enhances-Klebsiella-pneumoniae-pathogenesis
Source data and code for the associated manuscript

Hi there! Thanks for taking a look at this. This repository intends to store R files and raw data files associated with the above project. The goal of this repository is to provide open-access to the analysis performed for the above project, such that it can be replicated by any end-user. If you decide to run this analysis yourself, make sure to change all of the directories!

The data processing files are:

mothur_processing_dbdb.batch (data processing with mothur)

otu_corr.R (elastic net processing with mikropml using OTUs as input data)

combine_otu_corr_en.R (assembly of elastic net data using OTUs as input data)

Makefile (makefile rule structure for making .Rds files from mikropml)

Note For machine learning models in mikropml, I am providing a single data input example; however, I am happy to provide any additional set of processing files.

The analysis files are:

otu_dbdb.R (analysis and visualization of mothur and mikropml output)

asv_dbdb.R (analysis and visualization of mothur output)

amino_acid_analysis.R (analysis and visualization of metabolomics data)

in_vivo_analysis.R (analysis and visualization of mouse experimental data)

luminex.R (analysis and visualization of luminex data)

Finally, I am providing a mothur logfile (mothur.1715352818.logile) for the 16S rRNA gene sequencing data processing for user reference. 

