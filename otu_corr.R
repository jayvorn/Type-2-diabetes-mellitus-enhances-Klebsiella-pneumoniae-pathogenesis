#!/usr/bin/env Rscript

library(readxl)
library(tidyverse)
library(mikropml)
library(glue)
library(dplyr)

seed<-as.numeric(commandArgs(trailingOnly=TRUE))
#READ IN DATA##############
otu_metadata<-read.csv("raw_data/og_metadata.csv") %>%
  na.omit()

otu_counts<-read_tsv(file=glue("raw_data/final_dbdb.opti_mcc.0.03.subsample.shared")) %>%
  select(-label, -numOtus) %>%
  rename(sample = Group) %>%
  pivot_longer(-sample, names_to="otu", values_to = "count") %>%
  mutate(sample = substr(sample, 1, 4))

taxonomy<-read_tsv(file=glue("raw_data/final_dbdb.opti_mcc.0.03.cons.taxonomy")) %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),
           sep=";") %>%
  mutate(pretty_otu = str_replace(string=otu,
                                  pattern="otu",
                                  replacement = "OTU"),
         genus = str_replace(string=genus,
                             pattern="(.*)",
                             replacement="*\\1*"),
         genus = str_replace(string=genus,
                             pattern="\\*(.*)_unclassified\\*",
                             replacement="Unclassified<br>*\\1*"),
         taxon = glue("{pretty_otu} {genus}")) %>%
  select(otu, taxon)

otu_rel_abun_corr<-inner_join(otu_metadata, otu_counts, by="sample") %>%
  inner_join(., taxonomy, by = "otu") %>%
  filter(group == "het_No" | group == "dbdb_No") %>%
  group_by(sample) %>%
  mutate(rel_abund = 100*(count / sum(count))) %>%
  ungroup() %>%
  group_by(group, taxon) %>%
  mutate(x = mean(rel_abund)) %>%
  filter(x>=0.5) %>%
  select(-count, -x) %>%
  mutate(group = factor(group,
                        levels=c("het_No","dbdb_No","het_Yes","dbdb_Yes")))

otu_corr_data<-otu_rel_abun_corr%>%
  select(sample, taxon, rel_abund) %>%
  left_join(., otu_metadata, by = "sample") %>%
  select(taxon, rel_abund, log_CFU_per_g_total) %>%
  pivot_wider(., names_from = taxon, values_from = rel_abund)

#PROCESSING DATA####################
processed_data<-preprocess_data(otu_corr_data, outcome_colname = "log_CFU_per_g_total")$dat_transformed

#WRITING MODEL FORMULA###############
ml<-function(seed){
  results<-run_ml(processed_data, 
                method = "glmnet",
                outcome_colname = "log_CFU_per_g_total",
                find_feature_importance = TRUE,
                seed = seed)
  saveRDS(results, file=glue("processed_data/otu_corr_en_results_{seed}.Rds"))
  }
map(c(1:100), ml)

#WRITING RESULTS##########
saveRDS(results, file=glue("processed_data/otu_corr_en_results_{seed}.Rds"))
