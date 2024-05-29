#!/usr/bin/env Rscript

library(plyr)
library(tidyverse)
library(mikropml)

rds_files <- list.files(path = "~/Desktop/mothur_dbdb/analysis/dbdb_glmnet_corr/processed_data/",
                        pattern = "otu_corr_en_results_(\\d*).Rds",
                        full.names = TRUE)

taxon_data_results <- map(rds_files, readRDS)

taxon_data_results %>%
  map_dfr(pluck, "feature_importance") %>%
  write_tsv("processed_data/otu_corr_en_feature_importance.tsv")

taxon_data_results %>%
  map(pluck, "performance") %>%
  rbind.fill() %>%
  write_tsv("processed_data/otu_corr_en_performance.tsv")

taxon_data_results %>%
  map(pluck, "trained_model") %>%
  map_dfr(pluck, "results") %>%
  write_tsv("processed_data/otu_corr_en_trained_model.tsv")
