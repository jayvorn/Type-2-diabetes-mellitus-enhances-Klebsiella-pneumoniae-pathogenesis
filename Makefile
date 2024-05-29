#RULE STRUCTURE
#Rule
#target:prerequisites
#(tab)recipe

processed_data/otu_corr_en_results_%.Rds	:	raw_data/final_dbdb.opti_mcc.0.03.cons.taxonomy\
												raw_data//final_dbdb.opti_mcc.0.03.subsample.shared\
												raw_data/og_metadata.csv
	./code/otu_corr.R $*
	
SEEDS = $(shell seq 1 1 100)
EN_RDS = $(patsubst %,processed_data/otu_corr_en_results_%.Rds,$(SEEDS))

processed_data/otu_corr_en_%.tsv : code/combine_otu_corr_en.R $(EN_RDS)
	$^
