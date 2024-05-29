library(readxl)
library(tidyverse)
library(glue)
library(dplyr)
library(ggtext)
library(vegan)
library(ggplot2)
library(broom)
library(purrr)
library(mikropml)

#===================================================================
# SET ENVIRONMENT & IMPORT DATA
#===================================================================
# set path for data
setwd("~/Desktop/mothur_dbdb/")
# assign tag prefixes
dbdb_prefix<-"final_dbdb.opti_mcc"
# assign aesthetics
genotype_4_levels<-c("het_No","dbdb_No","het_Yes","dbdb_Yes")
genotype_4_labels<-c("db/+<br>No abx","db/db<br>No abx","db/+<br>Abx","db/db<br>Abx")
genotype_4_shapes<-c(21,24, 21,24)
genotype_4_colors<-c("black","black","#5F4B8B","#AD5E99")
genotype_4_fills<-c("#5F4B8B","#AD5E99","grey","lightgrey")
genotype_2_no_abx_levels<-c("het_No","dbdb_No")
genotype_2_no_abx_labels<-c("db/+<br>No abx","db/db<br>No abx")
genotype_2_no_abx_shapes<-c(21,24)
genotype_2_no_abx_colors<-c("black","black")
genotype_2_no_abx_fills<-c("#5F4B8B","#AD5E99")
genotype_2_abx_levels<-c("het_Yes","dbdb_Yes")
genotype_2_abx_labels<-c("db/+<br>Abx","db/db<br>Abx")
genotype_2_abx_shapes<-c(21,24)
genotype_2_abx_colors<-c("#5F4B8B","#AD5E99")
genotype_2_abx_fills<-c("grey","lightgrey")
# define sample size
run_samplesize<-function(tag){
  samplesize<-glue('./mothur "#count.groups(shared={tag}.shared, inputdir=raw_data)"')
  system(samplesize)
}
run_samplesize(dbdb_prefix)
# define subsample size
subsample_size<-3793
# create subsample function
run_subsample<-function(tag, size){
  subsample<-glue('./mothur "#sub.sample(shared={tag}.shared, size={size}, inputdir=raw_data, outputdir=raw_data)"')
  system(subsample)
}
# running subsample
run_subsample(dbdb_prefix, subsample_size)
# read and join otu and tax files
otu_metadata<-read.csv("~/Desktop/mothur_dbdb/raw_data/og_metadata.csv") %>%
  na.omit()
otu_counts<-read_tsv(file=glue("raw_data/{dbdb_prefix}.0.03.subsample.shared")) %>%
  select(-label, -numOtus) %>%
  rename(sample = Group) %>%
  pivot_longer(-sample, names_to="otu", values_to = "count") %>%
  mutate(sample = substr(sample, 1, 4))
taxonomy<-read_tsv(file=glue("raw_data/{dbdb_prefix}.0.03.cons.taxonomy")) %>%
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
complete_taxonomy<-read_tsv(file=glue("raw_data/{dbdb_prefix}.0.03.cons.taxonomy")) %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),
           sep=";")
otu_rel_abun<-inner_join(otu_metadata, otu_counts, by="sample") %>%
  inner_join(., taxonomy, by = "otu") %>%
  group_by(sample) %>%
  mutate(rel_abund = 100*(count / sum(count))) %>%
  ungroup() %>%
  select(-count) %>%
  mutate(group = factor(group,
                        levels=c("het_No","dbdb_No","het_Yes","dbdb_Yes")
  ))
otu_rel_abun_taxon_level<-inner_join(otu_metadata, otu_counts, by="sample") %>%
  inner_join(., complete_taxonomy, by = "otu") %>%
  group_by(sample) %>%
  mutate(rel_abund = 100*(count / sum(count))) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(c("kingdom", "phylum", "class", "order", "family", "genus", "otu"),
               names_to="level",
               values_to="taxon")
#===================================================================
# BETA DIVERSITY
#===================================================================
# generate data for dist matrix
dist_data_all<-otu_counts%>%
  pivot_wider(names_from = otu, values_from = count)%>%
  column_to_rownames("sample")
dist_data_noabx<-otu_counts%>%
  left_join(., otu_metadata, by = "sample") %>%
  filter(group == "het_No" | group == "dbdb_No") %>%
  select(sample, otu, count) %>%
  pivot_wider(names_from = otu, values_from = count)%>%
  column_to_rownames("sample")
dist_data_em08<-otu_counts%>%
  left_join(., otu_metadata, by = "sample") %>%
  filter(group == "het_No" | group == "dbdb_No") %>%
  filter(barrier == "em08") %>%
  select(sample, otu, count) %>%
  pivot_wider(names_from = otu, values_from = count)%>%
  column_to_rownames("sample")
dist_data_mp16<-otu_counts%>%
  left_join(., otu_metadata, by = "sample") %>%
  filter(group == "het_No" | group == "dbdb_No") %>%
  filter(barrier == "mp16") %>%
  select(sample, otu, count) %>%
  pivot_wider(names_from = otu, values_from = count)%>%
  column_to_rownames("sample")
dist_data_abx<-otu_counts%>%
  left_join(., otu_metadata, by = "sample") %>%
  filter(group == "het_Yes" | group == "dbdb_Yes") %>%
  select(sample, otu, count) %>%
  pivot_wider(names_from = otu, values_from = count)%>%
  column_to_rownames("sample")

# generate distance matrices
otu_dist_all=avgdist(dist_data_all, sample = subsample_size, dmethod = "robust.aitchison")
otu_dist_noabx=avgdist(dist_data_noabx, sample = subsample_size, dmethod = "robust.aitchison")
otu_dist_no_abx_em08=avgdist(dist_data_em08, sample = subsample_size, dmethod = "robust.aitchison")
otu_dist_no_abx_mp16=avgdist(dist_data_mp16, sample = subsample_size, dmethod = "robust.aitchison")
otu_dist_abx=avgdist(dist_data_abx, sample = subsample_size, dmethod = "robust.aitchison")

# generate amova metadata
amova_metadata_all<-otu_metadata %>%
  filter(sample != "M303" & sample != "M689" & sample != "M690" & sample != "M693")
amova_metadata_no_abx<-otu_metadata %>%
  filter(sample != "M303" & sample != "M689" & sample != "M690" & sample != "M693") %>%
  filter(group != "het_Yes" & group != "dbdb_Yes")
amova_metadata_no_abx_em08<-otu_metadata %>%
  filter(sample != "M303" & sample != "M689" & sample != "M690" & sample != "M693") %>%
  filter(group != "het_Yes" & group != "dbdb_Yes")%>%
  filter(barrier == "em08")
amova_metadata_no_abx_mp16<-otu_metadata %>%
  filter(sample != "M303" & sample != "M689" & sample != "M690" & sample != "M693") %>%
  filter(group != "het_Yes" & group != "dbdb_Yes")%>%
  filter(barrier == "mp16")
amova_metadata_abx<-otu_metadata %>%
  filter(sample != "M303" & sample != "M689" & sample != "M690" & sample != "M693")%>%
  filter(group != "het_No" & group != "dbdb_No")

# amova comparisons
aov_all<-adonis2(otu_dist_all ~ amova_metadata_all$group, permutations = 1000)
aov_no_abx<-adonis2(otu_dist_noabx ~ amova_metadata_no_abx$group, permutations = 1000)
aov_no_abx_em08<-adonis2(otu_dist_no_abx_em08 ~ amova_metadata_no_abx_em08$genotype, permutations = 1000)
aov_no_abx_mp16<-adonis2(otu_dist_no_abx_mp16 ~ amova_metadata_no_abx_mp16$genotype, permutations = 1000)
aov_abx<-adonis2(otu_dist_abx ~ amova_metadata_abx$genotype, permutations = 1000)

# plot pcoas
otu_all_pcoa=cmdscale(otu_dist_all, eig=TRUE, add=TRUE)
otu_all_pcoa_positions<-otu_all_pcoa$points
colnames(otu_all_pcoa_positions) = c("axis_1", "axis_2")

otu_all_pcoa_positions %>%
  as_tibble(rownames = "sample") %>%
  left_join(., otu_metadata, by = "sample") %>%
  select(group, sample, genotype, barrier, axis_1, axis_2) %>%
  ggplot(aes(x=axis_1, y=axis_2, color=group, shape = group, fill= group)) +
  geom_point(alpha = 1, size = 2.5) +
  stat_ellipse(level=0.7,
               geom="polygon",
               alpha=0.1,
               show.legend=FALSE) +
  labs(x = "Axis 1 (15.5%)", y = "Axis 2 (8.55%)") +
  scale_color_manual(name = "Group",
                     breaks = genotype_4_levels,
                     values = genotype_4_colors,
                     labels = genotype_4_labels) +
  scale_fill_manual(name = "Group",
                    breaks = genotype_4_levels,
                    values = genotype_4_fills,
                    labels = genotype_4_labels) +
  scale_shape_manual(name = "Group",
                     breaks = genotype_4_levels,
                     values = genotype_4_shapes,
                     labels = genotype_4_labels) +
  scale_y_continuous(limits = c(-13, 13),
                     breaks = c(-10,-5,0,5,10),
                     labels = c(-10,-5,0,5,10)) +
  scale_x_continuous(limits = c(-13, 13),
                     breaks = c(-10,-5,0,5,10),
                     labels = c(-10,-5,0,5,10)) +
  theme_classic() +
  theme(axis.title = element_markdown(size = 14),
        axis.text.x = element_markdown(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10))
ggsave("graphs/otu_PCoA_all.pdf", height = 15, width = 14, unit = "cm")

otu_no_abx_em08_pcoa=cmdscale(otu_dist_no_abx_em08, eig=TRUE, add=TRUE)
otu_no_abx_em08_pcoa_positions<-otu_no_abx_em08_pcoa$points
colnames(otu_no_abx_em08_pcoa_positions) = c("axis_1", "axis_2")

otu_no_abx_em08_pcoa_positions %>%
  as_tibble(rownames = "sample") %>%
  left_join(., amova_metadata_no_abx_em08, by = "sample") %>%
  select(group, sample, genotype, barrier, axis_1, axis_2) %>%
  ggplot(aes(x=axis_1, y=axis_2, color=group, shape = group, fill= group)) +
  geom_point(alpha = 1, size = 2.5) +
  stat_ellipse(level=0.6,
               geom="polygon",
               alpha=0.1,
               show.legend=FALSE) +
  labs(x = "Axis 1 (18.7%)", y = "Axis 2 (9.98%)") +
  scale_color_manual(name = "Group",
                     breaks = genotype_2_no_abx_levels,
                     values = genotype_2_no_abx_colors,
                     labels = genotype_2_no_abx_labels) +
  scale_fill_manual(name = "Group",
                    breaks = genotype_2_no_abx_levels,
                    values = genotype_2_no_abx_fills,
                    labels = genotype_2_no_abx_labels) +
  scale_shape_manual(name = "Group",
                     breaks = genotype_2_no_abx_levels,
                     values = genotype_2_no_abx_shapes,
                     labels = genotype_2_no_abx_labels) +
  scale_y_continuous(limits = c(-10, 10),
                     breaks = c(-10,-5,0,5,10),
                     labels = c(-10,-5,0,5,10)) +
  scale_x_continuous(limits = c(-15, 10),
                     breaks = c(-15,-10,-5,0,5,10),
                     labels = c(-15,-10,-5,0,5,10)) +
  theme_classic() +
  theme(axis.title = element_markdown(size = 14),
        axis.text.x = element_markdown(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10))
ggsave("graphs/otu_PCoA_no_abx_em08.pdf", height = 15, width = 14, unit = "cm")

otu_no_abx_mp16_pcoa=cmdscale(otu_dist_no_abx_mp16, eig=TRUE, add=TRUE)
otu_no_abx_mp16_pcoa_positions<-otu_no_abx_mp16_pcoa$points
colnames(otu_no_abx_mp16_pcoa_positions) = c("axis_1", "axis_2")

otu_no_abx_mp16_pcoa_positions %>%
  as_tibble(rownames = "sample") %>%
  left_join(., amova_metadata_no_abx_mp16, by = "sample") %>%
  select(group, sample, genotype, barrier, axis_1, axis_2) %>%
  ggplot(aes(x=axis_1, y=axis_2, color=group, shape = group, fill= group)) +
  geom_point(alpha = 1, size = 2.5) +
  stat_ellipse(level=0.6,
               geom="polygon",
               alpha=0.1,
               show.legend=FALSE) +
  labs(x = "Axis 1 (13.9%)", y = "Axis 2 (9.61%)") +
  scale_color_manual(name = "Group",
                     breaks = genotype_2_no_abx_levels,
                     values = genotype_2_no_abx_colors,
                     labels = genotype_2_no_abx_labels) +
  scale_fill_manual(name = "Group",
                    breaks = genotype_2_no_abx_levels,
                    values = genotype_2_no_abx_fills,
                    labels = genotype_2_no_abx_labels) +
  scale_shape_manual(name = "Group",
                     breaks = genotype_2_no_abx_levels,
                     values = genotype_2_no_abx_shapes,
                     labels = genotype_2_no_abx_labels) +
  scale_y_continuous(limits = c(-12, 12),
                     breaks = c(-10,-5,0,5,10),
                     labels = c(-10,-5,0,5,10)) +
  scale_x_continuous(limits = c(-10, 12),
                     breaks = c(-10,-5,0,5,10, 15),
                     labels = c(-10,-5,0,5,10, 15)) +
  theme_classic() +
  theme(axis.title = element_markdown(size = 14),
        axis.text.x = element_markdown(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10))
ggsave("graphs/otu_PCoA_no_abx_mp16.pdf", height = 15, width = 14, unit = "cm")

otu_abx_pcoa=cmdscale(otu_dist_abx, eig=TRUE, add=TRUE)
otu_abx_pcoa_positions<-otu_abx_pcoa$points
colnames(otu_abx_pcoa_positions) = c("axis_1", "axis_2")

otu_abx_pcoa_positions %>%
  as_tibble(rownames = "sample") %>%
  left_join(., amova_metadata_abx, by = "sample") %>%
  select(group, sample, genotype, barrier, axis_1, axis_2) %>%
  ggplot(aes(x=axis_1, y=axis_2, color=group, shape = group, fill= group)) +
  geom_point(alpha = 1, size = 2.5) +
  stat_ellipse(level=0.7,
               geom="polygon",
               alpha=0.1,
               show.legend=FALSE) +
  labs(x = "Axis 1 (35.2%)", y = "Axis 2 (19.8%)") +
  scale_color_manual(name = "Group",
                     breaks = genotype_2_abx_levels,
                     values = genotype_2_abx_colors,
                     labels = genotype_2_no_abx_labels) +
  scale_fill_manual(name = "Group",
                    breaks = genotype_2_abx_levels,
                    values = genotype_2_abx_fills,
                    labels = genotype_2_abx_labels) +
  scale_shape_manual(name = "Group",
                     breaks = genotype_2_abx_levels,
                     values = genotype_2_abx_shapes,
                     labels = genotype_2_abx_labels) +
  scale_y_continuous(limits = c(-10, 10),
                     breaks = c(-10,-5,0,5,10),
                     labels = c(-10,-5,0,5,10)) +
  scale_x_continuous(limits = c(-10, 20),
                     breaks = c(-10,-5,0,5,10,20),
                     labels = c(-10,-5,0,5,10,20)) +
  theme_classic() +
  theme(axis.title = element_markdown(size = 14),
        axis.text.x = element_markdown(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10))
ggsave("graphs/otu_PCoA_abx.pdf", height = 15, width = 14, unit = "cm")
#===================================================================
# LEFSE/LDA
#===================================================================
# No abx LefSe
shared_file<-read_tsv(file=glue("raw_data/{dbdb_prefix}.0.03.subsample.shared")) %>%
  rename(sample = Group)
# join share and metadata
shared_design<-inner_join(shared_file, otu_metadata, by="sample") %>%
  rename(microbiota = group) %>%
  rename(group = sample) %>%
# modify for barrier-specific LefSe  
  #filter(barrier == "mp16") %>%
  select(-X, -genotype, -barrier, -log_CFU_per_g_total)
# LefSe function
run_lefse<-function(x, y, tag){
  # select specific comparison groups
  x_y <-shared_design %>%
    filter(microbiota == x | microbiota == y)
  # write shared file for pairwise comparison
  x_y %>%
    select(-microbiota) %>%
    write_tsv(glue("processed_data/{tag}.shared"))
  # write design file for pairwise comparison
  x_y %>%
    select(group, microbiota) %>%
    write_tsv(glue("processed_data/{tag}.design"))
  # assign lefse function
  lefse<-glue('./mothur "#lefse(shared=processed_data/{tag}.shared, design=processed_data/{tag}.design, inputdir=processed_data)"')
  # Run lefse function
  system(lefse)
  return(glue("processed_data/{tag}.0.03.lefse_summary"))
}
# run LefSe
NoAbx<-run_lefse("het_No","dbdb_No","noabx")
Barrier_em08<-run_lefse("het_No","dbdb_No","em08")
Barrier_mp16<-run_lefse("het_No","dbdb_No","mp16")
# LDA Graphs
# noabx lefse
read_tsv(NoAbx) %>%
  drop_na(LDA) %>%
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  filter(LDA >= 3.5) %>%
  mutate(LDA = if_else(Class == "het_No", -1*LDA, LDA),
         taxon = fct_reorder(taxon, LDA),
         label_x = if_else(Class == "het_No", 0.1, -0.1),
         label_hjust = if_else(Class == "het_No", 0, 1)) %>%
  ggplot(aes(x=LDA, y=taxon, label=taxon, fill=Class)) +
  geom_col() +
  geom_richtext(aes(x=label_x, hjust=label_hjust), fill = NA, label.color = NA, size = 2.5) + 
  labs(y=NULL, x=bquote("log"[10]*"(LDA Score)")) +
  scale_fill_manual(name = "Group",
                    breaks = genotype_2_no_abx_levels,
                    values = genotype_2_no_abx_fills,
                    labels = genotype_2_no_abx_labels) +
  scale_x_continuous(limits=c(-5,5))+
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown())
ggsave(glue("graphs/noabx_otu_LefSe.pdf"), height = 2.5, width = 3.5)
# em08 lefse
read_tsv(Barrier_em08) %>%
  drop_na(LDA) %>%
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  filter(LDA >= 3.5) %>%
  mutate(LDA = if_else(Class == "het_No", -1*LDA, LDA),
         taxon = fct_reorder(taxon, LDA),
         label_x = if_else(Class == "het_No", 0.1, -0.1),
         label_hjust = if_else(Class == "het_No", 0, 1)) %>%
  ggplot(aes(x=LDA, y=taxon, label=taxon, fill=Class)) +
  geom_col() +
  labs(title = "Barrier em08") +
  geom_richtext(aes(x=label_x, hjust=label_hjust), fill = NA, label.color = NA, size = 2.5) + 
  labs(y=NULL, x=bquote("log"[10]*"(LDA Score)")) +
  scale_fill_manual(name = "Group",
                    breaks = genotype_2_no_abx_levels,
                    values = genotype_2_no_abx_fills,
                    labels = genotype_2_no_abx_labels) +
  scale_x_continuous(limits=c(-5.5,5.5))+
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown())
ggsave(glue("graphs/barrier_em08_otu_LefSe.pdf"), height = 4, width = 4)
# mp16 lefse
read_tsv(Barrier_mp16) %>%
  drop_na(LDA) %>%
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  filter(LDA >= 3.5) %>%
  mutate(LDA = if_else(Class == "het_No", -1*LDA, LDA),
         taxon = fct_reorder(taxon, LDA),
         label_x = if_else(Class == "het_No", 0.1, -0.1),
         label_hjust = if_else(Class == "het_No", 0, 1)) %>%
  ggplot(aes(x=LDA, y=taxon, label=taxon, fill=Class)) +
  geom_col() +
  labs(title = "Barrier mp16") +
  geom_richtext(aes(x=label_x, hjust=label_hjust), fill = NA, label.color = NA, size = 2.5) + 
  labs(y=NULL, x=bquote("log"[10]*"(LDA Score)")) +
  scale_fill_manual(name = "Group",
                    breaks = genotype_2_no_abx_levels,
                    values = genotype_2_no_abx_fills,
                    labels = genotype_2_no_abx_labels) +
  scale_x_continuous(limits=c(-5.5,5.5))+
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown())
ggsave(glue("graphs/barrier_mp16_otu_LefSe.pdf"), height = 4, width = 4)
#===================================================================
# ALPHA DIVERSITY
#===================================================================
#Write alpha diversity function
run_alpha<-function(tag){
  alpha<-glue('./mothur "#summary.single(shared={tag}.0.03.subsample.shared, subsample=F, inputdir=raw_data, outputdir=raw_data)"')
  system(alpha)
}

#Running alpha diversity function
run_alpha(dbdb_prefix)

#Create alpha diversity metadata file
alpha_diversity <- read_tsv(glue("raw_data/{dbdb_prefix}.0.03.subsample.groups.summary")) %>%
  mutate(invsimpson = 1/simpson)%>%
  mutate(group = substr(group, 1, 4)) %>%
  rename(sample = group)
metadata_alpha <- inner_join(otu_metadata, alpha_diversity, by="sample")

#Test for differences in alpha diversity
sha_aov<-TukeyHSD(aov(log(shannon) ~ group, data = metadata_alpha))
chao_aov<-TukeyHSD(aov(log(chao) ~ group, data = metadata_alpha))
IS_aov<-TukeyHSD(aov(log(invsimpson) ~ group, data = metadata_alpha))
sobs_aov<-TukeyHSD(aov(log(sobs) ~ group, data = metadata_alpha))
#Graph inverse simpson
metadata_alpha %>%
  ggplot(aes(x=group, y=invsimpson, shape=group, color=group, fill=group)) +
  geom_boxplot(show.legend = FALSE, width = 0.5, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             show.legend = FALSE) +
  labs(x=NULL, 
       y="Inverse Simpson") +
  scale_color_manual(name = "Group",
                     breaks = genotype_4_levels,
                     values = genotype_4_colors,
                     labels = genotype_4_labels) +
  scale_fill_manual(name = "Group",
                    breaks = genotype_4_levels,
                    values = genotype_4_fills,
                    labels = genotype_4_labels) +
  scale_shape_manual(name = "Group",
                     breaks = genotype_4_levels,
                     values = genotype_4_shapes,
                     labels = genotype_4_labels) +
  scale_y_continuous(limits = c(0, 46),
                     breaks = c(0,10,20,30,40),
                     labels = c(0,10,20,30,40))+
  scale_x_discrete(breaks = genotype_4_levels,
                   labels = genotype_4_labels)+
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 12),
        legend.title = element_blank()) +
  annotate("segment", x = c(1,1,3), xend = c(2,3,4), 
           y = c(37,43,37), yend = c(37,43,37), 
           color = "black")+
  geom_richtext(data=tibble(x=1.5, y=40), fill = NA, label.color = NA, label="*p* < 1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=2, y=46), fill = NA, label.color = NA, label="*p* = 0.33" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=3.5, y=40), fill = NA, label.color = NA, label="*p* < 1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=3)
ggsave(glue("graphs/dbdb_otu_inverse_simpson.pdf"), height = 6, width = 8, unit = "cm")

#Graph shannon
metadata_alpha %>%
  ggplot(aes(x=group, y=shannon, shape=group, color=group, fill=group)) +
  geom_boxplot(show.legend = FALSE, width = 0.5, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             show.legend = FALSE) +
  labs(x=NULL, 
       y="Shannon") +
  scale_color_manual(name = "Group",
                     breaks = genotype_4_levels,
                     values = genotype_4_colors,
                     labels = genotype_4_labels) +
  scale_fill_manual(name = "Group",
                    breaks = genotype_4_levels,
                    values = genotype_4_fills,
                    labels = genotype_4_labels) +
  scale_shape_manual(name = "Group",
                     breaks = genotype_4_levels,
                     values = genotype_4_shapes,
                     labels = genotype_4_labels) +
  scale_y_continuous(limits = c(0, 5.575),
                     breaks = c(0,1,2,3,4,5),
                     labels = c(0,1,2,3,4,5))+
  scale_x_discrete(breaks = genotype_4_levels,
                   labels = genotype_4_labels)+
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 12),
        legend.title = element_blank()) +
  annotate("segment", x = c(1,1,3), xend = c(2,3,4), 
           y = c(4.6,5.3,4.6), yend = c(4.6,5.3,4.6), 
           color = "black")+
  geom_richtext(data=tibble(x=1.5, y=4.875), fill = NA, label.color = NA, label="*p* = 6.0e-4" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=2, y=5.575), fill = NA, label.color = NA, label="*p* = 0.99" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=3.5, y=4.875), fill = NA, label.color = NA, label="*p* < 1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=3)
ggsave(glue("graphs/dbdb_otu_shannon.pdf"), height = 6, width = 8, unit = "cm")

#Graph chao
metadata_alpha %>%
  ggplot(aes(x=group, y=chao, shape=group, color=group, fill=group)) +
  geom_boxplot(show.legend = FALSE, width = 0.5, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             show.legend = FALSE) +
  labs(x=NULL, 
       y="Chao") +
  scale_color_manual(name = "Group",
                     breaks = genotype_4_levels,
                     values = genotype_4_colors,
                     labels = genotype_4_labels) +
  scale_fill_manual(name = "Group",
                    breaks = genotype_4_levels,
                    values = genotype_4_fills,
                    labels = genotype_4_labels) +
  scale_shape_manual(name = "Group",
                     breaks = genotype_4_levels,
                     values = genotype_4_shapes,
                     labels = genotype_4_labels) +
  scale_y_continuous(limits = c(0, 480),
                     breaks = c(0,100,200,300,400),
                     labels = c(0,100,200,300,400))+
  scale_x_discrete(breaks = genotype_4_levels,
                   labels = genotype_4_labels)+
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 12),
        legend.title = element_blank()) +
  annotate("segment", x = c(1,1,3), xend = c(2,3,4), 
           y = c(350,450,400), yend = c(350,450,400), 
           color = "black")+
  geom_richtext(data=tibble(x=1.5, y=380), fill = NA, label.color = NA, label="*p* < 1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=2, y=480), fill = NA, label.color = NA, label="*p* = 0.99" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=3.5, y=430), fill = NA, label.color = NA, label="*p* < 1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=3)
ggsave(glue("graphs/dbdb_otu_chao.pdf"), height = 6, width = 8, unit = "cm")
#===================================================================
# CORR MIKROPML
#===================================================================
# create input data
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

# process data
processed_data<-preprocess_data(otu_corr_data, outcome_colname = "log_CFU_per_g_total")$dat_transformed

# model forumula
ml<-function(seed){
  results<-run_ml(processed_data, 
                  method = "glmnet",
                  outcome_colname = "log_CFU_per_g_total",
                  find_feature_importance = TRUE,
                  seed = seed)
  saveRDS(results, file=glue("processed_data/otu_corr_en_results_{seed}.Rds"))
}
map(c(1:100), ml)

#Create en file list
rds_files <- list.files(path = "~/Desktop/mothur_dbdb/analysis/dbdb_glmnet_corr/processed_data/",
                        pattern = "otu_corr_en_results_(\\d*).Rds",
                        full.names = TRUE)

taxon_data_results <- map(rds_files, readRDS)

taxon_data_results %>%
  map_dfr(pluck, "feature_importance") %>%
  write_tsv("~/Desktop/mothur_dbdb/analysis/dbdb_glmnet_corr/processed_data/otu_corr_en_feature_importance.tsv")

#library(plyr)
#remember to unload plyr before moving on
taxon_data_results %>%
  map(pluck, "performance") %>%
  rbind.fill() %>%
  write_tsv("~/Desktop/mothur_dbdb/analysis/dbdb_glmnet_corr/processed_data/otu_corr_en_performance.tsv")

taxon_data_results %>%
  map(pluck, "trained_model") %>%
  map_dfr(pluck, "results") %>%
  write_tsv("~/Desktop/mothur_dbdb/analysis/dbdb_glmnet_corr/processed_data/otu_corr_en_trained_model.tsv")

#Write function to extract en asv feature weights
get_weights<-function(file_name){
  model<-readRDS(file_name) %>%
    pluck("trained_model")
  coef(model$finalModel, model$bestTune$lambda) %>%
    as.matrix() %>%
    as.tibble(rownames = "feature") %>%
    rename(weight = s1) %>%
    mutate(seed = str_replace(file_name, "~/Desktop/mothur_dbdb/analysis/dbdb_glmnet_corr/processed_data/otu_corr_en_results_(\\d*).Rds", "\\1"))
}

model<-readRDS("/Users/jayvornh/Desktop/mothur_dbdb/analysis/dbdb_glmnet_corr/processed_data//otu_corr_en_results_96.Rds")%>%
  pluck("trained_model")

X<-coef(model$finalModel, model$bestTune$lambda)%>%
  as.matrix() %>%
  as.tibble(rownames = "feature") %>%
  rename(weight = s1) %>%

#Run en feature weight extraction function
weights<-map_dfr(rds_files, get_weights)

#Remove intercept variable and categorize barrier status
weight_summary<-weights %>%
  group_by(feature) %>%
  summarize(mean_weight = mean(weight)) %>%
  filter(feature != "(Intercept)") %>%
  mutate(feature = str_remove(feature, "^`")) %>%
  mutate(feature = str_remove(feature, "`$"))

feature_importance<-read_tsv(glue("~/Desktop/mothur_dbdb/analysis/dbdb_glmnet_corr/processed_data/otu_corr_en_feature_importance.tsv"))
feature_importance %>%
  rename(feature = feat) %>%
  group_by(feature) %>%
  summarize(med_perf_metric_diff = median(perf_metric_diff * 1),
            l_quartile_perf_metric_diff = quantile((perf_metric_diff) * 1, prob = 0.25),
            u_quartile_perf_metric_diff = quantile((perf_metric_diff) * 1, prob = 0.75)) %>%
  inner_join(., weight_summary, by = "feature") %>%
  filter(med_perf_metric_diff < -0.025) %>%
  mutate(feature = fct_reorder(feature, med_perf_metric_diff)) %>%
  ggplot(aes(x=med_perf_metric_diff/-1, y=feature, xmin=l_quartile_perf_metric_diff/-1, xmax = u_quartile_perf_metric_diff/-1)) +
  geom_vline(xintercept=0, linetype = "dashed", color = "darkgrey") +
  geom_point(show.legend = FALSE, size = 3) +
  geom_linerange(show.legend = FALSE) +
  xlim(-0.01, 0.08) +
  labs(x = "Increase in RSME") +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 8),
        axis.title = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        axis.text = element_markdown(color = "black", size = 10),
        axis.text.y = element_markdown(),
        axis.title.y = element_blank(),
        legend.title = element_blank())
ggsave(glue("graphs/dbdb_rsme_diff.pdf"), height = 5.75, width = 10, unit = "cm")

#build r plot
corr_data_Otu0006<-otu_rel_abun %>%
  filter(otu == "Otu0006") %>%
  select(log_CFU_per_g_total, rel_abund)
Otu0006_r<-cor(corr_data_Otu0006$rel_abund, corr_data_Otu0006$log_CFU_per_g_total, method = "spearman")

corr_data_Otu0010<-otu_rel_abun %>%
  filter(otu == "Otu0021") %>%
  select(log_CFU_per_g_total, rel_abund)
Otu0010_r<-cor(corr_data_Otu0010$rel_abund, corr_data_Otu0010$log_CFU_per_g_total, method = "spearman")

corr_data_Otu0013<-otu_rel_abun %>%
  filter(otu == "Otu0013") %>%
  select(log_CFU_per_g_total, rel_abund)
Otu0013_r<-cor(corr_data_Otu0013$rel_abund, corr_data_Otu0013$log_CFU_per_g_total, method = "spearman")

corr_data_Otu0015<-otu_rel_abun %>%
  filter(otu == "Otu0021") %>%
  select(log_CFU_per_g_total, rel_abund)
Otu0015_r<-cor(corr_data_Otu0015$rel_abund, corr_data_Otu0015$log_CFU_per_g_total, method = "spearman")

corr_data_Otu0017<-otu_rel_abun %>%
  filter(otu == "Otu0017") %>%
  select(log_CFU_per_g_total, rel_abund)
Otu0017_r<-cor(corr_data_Otu0017$rel_abund, corr_data_Otu0017$log_CFU_per_g_total, method = "spearman")

corr_data_Otu0063<-otu_rel_abun %>%
  filter(otu == "Otu0063") %>%
  select(log_CFU_per_g_total, rel_abund)
Otu0063_r<-cor(corr_data_Otu0063$rel_abund, corr_data_Otu0063$log_CFU_per_g_total, method = "spearman")

as_tibble(rbind(Otu0006_r, Otu0013_r, Otu0017_r, Otu0063_r, Otu0015_r, Otu0010_r), rownames = "otu") %>%
  rename(r = V1) %>%
  mutate(otu = substr(otu, start = 1, stop = 7)) %>%
  left_join(otu_rel_abun,., by = "otu") %>%
  na.omit() %>%
  ggplot(aes(x = r, y = taxon, color = r)) +
  geom_vline(xintercept=0, linetype = "dashed", color = "darkgrey") +
  geom_point(show.legend = FALSE, size = 2.5) +
  scale_y_discrete(limits = c("Otu0006 *Paramuribaculum*", "Otu0013 Unclassified<br>*Muribaculaceae*", 
                              "Otu0017 *Dubosiella*", "Otu0063 Unclassified<br>*Lachnospiraceae*",
                              "Otu0015 Unclassified<br>*Muribaculaceae*", "Otu0010 Unclassified<br>*Muribaculaceae*"),
                   labels = c("Otu0006 *Paramuribaculum*", "Otu0013 Unclassified<br>*Muribaculaceae*", 
                              "Otu0017 *Dubosiella*", "Otu0063 *Clostridium_XIVa*",
                              "Otu0015 Unclassified<br>*Muribaculaceae*", "Otu0010 Unclassified<br>*Muribaculaceae*")) +
  scale_color_gradient(low = "#AD5E99", high = "#5F4B8B") +
  labs(x = "r") +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 8),
        axis.title = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        axis.text = element_markdown(color = "black", size = 10),
        axis.text.y = element_markdown(),
        axis.title.y = element_blank(),
        legend.title = element_blank())
ggsave(glue("graphs/dbdb_r_diff.pdf"), height = 5.75, width = 9, unit = "cm")
#===================================================================
# SIG ASVs
#===================================================================
#Determine significantly different ASVs
#Graphing otus
top_otus<-otu_rel_abun %>%
  filter(otu == "Otu0002" | otu == "Otu0006" | otu == "Otu0009" | otu == "Otu0010" | otu == "Otu0013" | 
           otu == "Otu0015" |otu == "Otu0017" |otu == "Otu0022" | otu == "Otu0034"| otu == "Otu0063")
top_otu_aov<-top_otus  %>%
  filter(otu == "Otu0063")
TukeyHSD(aov(rel_abund ~ group, data = top_otu_aov))

ggplot(top_otus, aes(y=otu, x=rel_abund, shape=group, color=group, fill=group)) +
  geom_boxplot(width = 0.8, size=0.3, alpha=0.25, outlier.shape = NA, coef= 1.5, fatten = 0.75) + 
  geom_point(position = position_dodge(width = 0.8), size = 2, alpha=1) + 
  labs(y=NULL, 
       x="Relative abundance (%)") +
  xlim(0, 55) +
  scale_color_manual(name = "Group",
                     breaks = genotype_4_levels,
                     values = genotype_4_colors,
                     labels = genotype_4_labels) +
  scale_fill_manual(name = "Group",
                    breaks = genotype_4_levels,
                    values = genotype_4_fills,
                    labels = genotype_4_labels) +
  scale_shape_manual(name = "Group",
                     breaks = genotype_4_levels,
                     values = genotype_4_shapes,
                     labels = genotype_4_labels) +
  theme_classic() +
  theme(legend.key.height = unit(0.75, "cm"),
        legend.position = "right",
        legend.margin = margin(t=1, b=2, r=2, l=2),
        legend.title = element_markdown(size=10),
        legend.text = element_markdown(size=6),
        axis.text.x = element_markdown(color = "black", size = 10),
        axis.text.y = element_markdown(color = "black",size = 10, angle = 270, hjust = 0.5)) +
  geom_richtext(data=tibble(x=52, y=0.8), fill = NA, label.color = NA, label="*p* = 0.03" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=25, y=1.8), fill = NA, label.color = NA, label="*p* = 0.94" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=20, y=2.8), fill = NA, label.color = NA, label="*p* = 0.045" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=17, y=3.8), fill = NA, label.color = NA, label="*p* = 3.2e-5" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=15, y=4.8), fill = NA, label.color = NA, label="*p* = 0.91" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=15, y=5.8), fill = NA, label.color = NA, label="*p* = 0.96" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=17, y=6.8), fill = NA, label.color = NA, label="*p* = 0.03" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=10, y=7.8), fill = NA, label.color = NA, label="*p* = 0.45" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=14, y=8.8), fill = NA, label.color = NA, label="*p* = 7.1e-4" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=10, y=9.8), fill = NA, label.color = NA, label="*p* = 0.84" ,aes(x=x, y=y), inherit.aes=FALSE, size=3)
ggsave("graphs/sig_otus.pdf", height = 8, width = 3.5)