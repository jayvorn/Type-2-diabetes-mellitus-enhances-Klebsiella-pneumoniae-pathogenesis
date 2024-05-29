library(readxl)
library(tidyverse)
library(glue)
library(dplyr)
library(ggtext)
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(broom)
library(purrr)

#===================================================================
# SET ENVIRONMENT & IMPORT DATA
#===================================================================
# set path for data
setwd("~/Desktop/mothur_dbdb/")
# assign tag prefixes
dbdb_prefix<-"final_dbdb.ASV"
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
# read and join asv and tax files
asv_metadata<-read.csv("~/Desktop/mothur_dbdb/raw_data/og_metadata.csv") %>%
  na.omit()
asv_counts<-read_tsv(file=glue("raw_data/{dbdb_prefix}.ASV.subsample.shared")) %>%
  select(-label, -numOtus) %>%
  rename(sample = Group) %>%
  pivot_longer(-sample, names_to="otu", values_to = "count") %>%
  mutate(sample = substr(sample, 1, 4))
taxonomy<-read_tsv(file=glue("raw_data/{dbdb_prefix}.ASV.cons.taxonomy")) %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),
           sep=";") %>%
  mutate(pretty_asv = str_replace(string=otu,
                                  pattern="ASV",
                                  replacement = "ASV"),
         genus = str_replace(string=genus,
                             pattern="(.*)",
                             replacement="*\\1*"),
         genus = str_replace(string=genus,
                             pattern="\\*(.*)_unclassified\\*",
                             replacement="Unclassified<br>*\\1*"),
         taxon = glue("{pretty_asv} {genus}")) %>%
  select(otu, taxon)
complete_taxonomy<-read_tsv(file=glue("raw_data/{dbdb_prefix}.ASV.cons.taxonomy")) %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),
           sep=";")
asv_rel_abun<-inner_join(asv_metadata, asv_counts, by="sample") %>%
  inner_join(., taxonomy, by = "otu") %>%
  group_by(group) %>%
  mutate(rel_abund = 100*(count / sum(count))) %>%
  ungroup() %>%
  select(-count) %>%
  mutate(group = factor(group,
                             levels=c("het _ No","dbdb _ No","het _ Yes","dbdb _ Yes")
  ))
asv_rel_abun_taxon_level<-inner_join(asv_metadata, asv_counts, by="sample") %>%
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
dist_data_all<-asv_counts%>%
  pivot_wider(names_from = otu, values_from = count)%>%
  column_to_rownames("sample")
dist_data_em08<-asv_counts%>%
  left_join(., asv_metadata, by = "sample") %>%
  filter(group == "het_No" | group == "dbdb_No") %>%
  filter(barrier == "em08") %>%
  select(sample, otu, count) %>%
  pivot_wider(names_from = otu, values_from = count)%>%
  column_to_rownames("sample")
dist_data_mp16<-asv_counts%>%
  left_join(., asv_metadata, by = "sample") %>%
  filter(group == "het_No" | group == "dbdb_No") %>%
  filter(barrier == "mp16") %>%
  select(sample, otu, count) %>%
  pivot_wider(names_from = otu, values_from = count)%>%
  column_to_rownames("sample")
dist_data_abx<-asv_counts%>%
  left_join(., asv_metadata, by = "sample") %>%
  filter(group == "het_Yes" | group == "dbdb_Yes") %>%
  select(sample, otu, count) %>%
  pivot_wider(names_from = otu, values_from = count)%>%
  column_to_rownames("sample")

# generate distance matrices
asv_dist_all=avgdist(dist_data_all, sample = subsample_size, dmethod = "robust.aitchison")
asv_dist_no_abx_em08=avgdist(dist_data_em08, sample = subsample_size, dmethod = "robust.aitchison")
asv_dist_no_abx_mp16=avgdist(dist_data_mp16, sample = subsample_size, dmethod = "robust.aitchison")
asv_dist_abx=avgdist(dist_data_abx, sample = subsample_size, dmethod = "robust.aitchison")

# generate amova metadata
amova_metadata_all<-asv_metadata %>%
  filter(sample != "M303" & sample != "M689" & sample != "M690" & sample != "M693")
amova_metadata_no_abx_em08<-asv_metadata %>%
  filter(sample != "M303" & sample != "M689" & sample != "M690" & sample != "M693") %>%
  filter(group != "het_Yes" & group != "dbdb_Yes")%>%
  filter(barrier == "em08")
amova_metadata_no_abx_mp16<-asv_metadata %>%
  filter(sample != "M303" & sample != "M689" & sample != "M690" & sample != "M693") %>%
  filter(group != "het_Yes" & group != "dbdb_Yes")%>%
  filter(barrier == "mp16")
amova_metadata_abx<-asv_metadata %>%
  filter(sample != "M303" & sample != "M689" & sample != "M690" & sample != "M693")%>%
  filter(group != "het_No" & group != "dbdb_No")

# amova comparisons
aov_all<-adonis2(asv_dist_all ~ amova_metadata_all$group, permutations = 1000)
aov_no_abx_em08<-adonis2(asv_dist_no_abx_em08 ~ amova_metadata_no_abx_em08$genotype, permutations = 1000)
aov_no_abx_mp16<-adonis2(asv_dist_no_abx_mp16 ~ amova_metadata_no_abx_mp16$genotype, permutations = 1000)
aov_abx<-adonis2(asv_dist_abx ~ amova_metadata_abx$genotype, permutations = 1000)

# plot pcoas
asv_all_pcoa=cmdscale(asv_dist_all, eig=TRUE, add=TRUE)
asv_all_pcoa_positions<-asv_all_pcoa$points
colnames(asv_all_pcoa_positions) = c("axis_1", "axis_2")

asv_all_pcoa_positions %>%
  as_tibble(rownames = "sample") %>%
  left_join(., asv_metadata, by = "sample") %>%
  select(group, sample, genotype, barrier, axis_1, axis_2) %>%
  ggplot(aes(x=axis_1, y=axis_2, color=group, shape = group, fill= group)) +
  geom_point(alpha = 1, size = 2.5) +
  stat_ellipse(level=0.8,
               geom="polygon",
               alpha=0.1,
               show.legend=FALSE) +
  labs(x = "Axis 1 (14.3%)", y = "Axis 2 (7.77%)") +
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
  scale_y_continuous(limits = c(-15, 15),
                     breaks = c(-15,-10,-5,0,5,10,15),
                     labels = c(-15,-10,-5,0,5,10,15)) +
  scale_x_continuous(limits = c(-15, 15),
                     breaks = c(-15,-10,-5,0,5,10,15),
                     labels = c(-15,-10,-5,0,5,10,15)) +
  theme_classic() +
  theme(axis.title = element_markdown(size = 14),
        axis.text.x = element_markdown(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10))
ggsave("graphs/asv_PCoA_all.pdf", height = 15, width = 14, unit = "cm")

asv_no_abx_em08_pcoa=cmdscale(asv_dist_no_abx_em08, eig=TRUE, add=TRUE)
asv_no_abx_em08_pcoa_positions<-asv_no_abx_em08_pcoa$points
colnames(asv_no_abx_em08_pcoa_positions) = c("axis_1", "axis_2")

asv_no_abx_em08_pcoa_positions %>%
  as_tibble(rownames = "sample") %>%
  left_join(., amova_metadata_no_abx_em08, by = "sample") %>%
  select(group, sample, genotype, barrier, axis_1, axis_2) %>%
  ggplot(aes(x=axis_1, y=axis_2, color=group, shape = group, fill= group)) +
  geom_point(alpha = 1, size = 2.5) +
  stat_ellipse(level=0.8,
               geom="polygon",
               alpha=0.1,
               show.legend=FALSE) +
  labs(x = "Axis 1 (16.3%)", y = "Axis 2 (11.1%)") +
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
  scale_y_continuous(limits = c(-15, 15),
                     breaks = c(-15,-10,-5,0,5,10, 15),
                     labels = c(-15,-10,-5,0,5,10, 15)) +
  scale_x_continuous(limits = c(-10, 17),
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
ggsave("graphs/asv_PCoA_no_abx_em08.pdf", height = 15, width = 14, unit = "cm")

asv_no_abx_mp16_pcoa=cmdscale(asv_dist_no_abx_mp16, eig=TRUE, add=TRUE)
asv_no_abx_mp16_pcoa_positions<-asv_no_abx_mp16_pcoa$points
colnames(asv_no_abx_mp16_pcoa_positions) = c("axis_1", "axis_2")

asv_no_abx_mp16_pcoa_positions %>%
  as_tibble(rownames = "sample") %>%
  left_join(., amova_metadata_no_abx_mp16, by = "sample") %>%
  select(group, sample, genotype, barrier, axis_1, axis_2) %>%
  ggplot(aes(x=axis_1, y=axis_2, color=group, shape = group, fill= group)) +
  geom_point(alpha = 1, size = 2.5) +
  stat_ellipse(level=0.6,
               geom="polygon",
               alpha=0.1,
               show.legend=FALSE) +
  labs(x = "Axis 1 (16.3%)", y = "Axis 2 (11.1%)") +
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
  scale_y_continuous(limits = c(-15, 15),
                     breaks = c(-15,-10,-5,0,5,10, 15),
                     labels = c(-15,-10,-5,0,5,10, 15)) +
  scale_x_continuous(limits = c(-10, 17),
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
ggsave("graphs/asv_PCoA_no_abx_mp16.pdf", height = 15, width = 14, unit = "cm")

asv_abx_pcoa=cmdscale(asv_dist_abx, eig=TRUE, add=TRUE)
asv_abx_pcoa_positions<-asv_abx_pcoa$points
colnames(asv_abx_pcoa_positions) = c("axis_1", "axis_2")

asv_abx_pcoa_positions %>%
  as_tibble(rownames = "sample") %>%
  left_join(., amova_metadata_abx, by = "sample") %>%
  select(group, sample, genotype, barrier, axis_1, axis_2) %>%
  ggplot(aes(x=axis_1, y=axis_2, color=group, shape = group, fill= group)) +
  geom_point(alpha = 1, size = 2.5) +
  stat_ellipse(level=0.7,
               geom="polygon",
               alpha=0.1,
               show.legend=FALSE) +
  labs(x = "Axis 1 (39.9%)", y = "Axis 2 (22.1%)") +
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
ggsave("graphs/asv_PCoA_abx.pdf", height = 15, width = 14, unit = "cm")
#===================================================================
# LEFSE/LDA
#===================================================================
# No abx LefSe
shared_file<-read_tsv(file=glue("raw_data/{dbdb_prefix}.ASV.subsample.shared")) %>%
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
  return(glue("processed_data/{tag}.ASV.lefse_summary"))
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
ggsave(glue("graphs/noabx_asv_LefSe.pdf"), height = 2.3, width = 3)
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
ggsave(glue("graphs/barrier_em08_asv_LefSe.pdf"), height = 4, width = 4)
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
ggsave(glue("graphs/barrier_mp16_asv_LefSe.pdf"), height = 4, width = 4)
#===================================================================
# ALPHA DIVERSITY
#===================================================================
#Write alpha diversity function
run_alpha<-function(tag){
  alpha<-glue('./mothur "#summary.single(shared={tag}.ASV.subsample.shared, subsample=F, inputdir=raw_data, outputdir=raw_data)"')
  system(alpha)
}

#Running alpha diversity function
run_alpha(dbdb_prefix)

#Create alpha diversity metadata file
alpha_diversity <- read_tsv(glue("raw_data/{dbdb_prefix}.ASV.subsample.groups.summary")) %>%
  mutate(invsimpson = 1/simpson)%>%
  mutate(group = substr(group, 1, 4)) %>%
  rename(sample = group)
metadata_alpha <- inner_join(asv_metadata, alpha_diversity, by="sample")

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
  scale_y_continuous(limits = c(0, 51),
                     breaks = c(0,10,20,30,40,50),
                     labels = c(0,10,20,30,40,50))+
  scale_x_discrete(breaks = genotype_4_levels,
                   labels = genotype_4_labels)+
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 12),
        legend.title = element_blank()) +
  annotate("segment", x = c(1,1,3), xend = c(2,3,4), 
           y = c(42,48,42), yend = c(42,48,42), 
           color = "black")+
  geom_richtext(data=tibble(x=1.5, y=44.5), fill = NA, label.color = NA, label="*p* < 1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=2, y=50.5), fill = NA, label.color = NA, label="*p* = 0.35" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=3.5, y=44.5), fill = NA, label.color = NA, label="*p* < 1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=3)
ggsave(glue("graphs/dbdb_asv_inverse_simpson.pdf"), height = 6, width = 8, unit = "cm")

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
  scale_y_continuous(limits = c(0, 5.8),
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
           y = c(4.7,5.5,4.7), yend = c(4.7,5.5,4.7), 
           color = "black")+
  geom_richtext(data=tibble(x=1.5, y=5), fill = NA, label.color = NA, label="*p* = 5.8e-4" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=2, y=5.8), fill = NA, label.color = NA, label="*p* = 0.99" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=3.5, y=5), fill = NA, label.color = NA, label="*p* < 1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=3)
ggsave(glue("graphs/dbdb_asv_shannon.pdf"), height = 6, width = 8, unit = "cm")

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
  scale_y_continuous(limits = c(0, 1000),
                     breaks = c(0,250,500,750,1000),
                     labels = c(0,250,500,750,1000))+
  scale_x_discrete(breaks = genotype_4_levels,
                   labels = genotype_4_labels)+
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 12),
        legend.title = element_blank()) +
  annotate("segment", x = c(1,1,3), xend = c(2,3,4), 
           y = c(800,950,800), yend = c(800,950,800), 
           color = "black")+
  geom_richtext(data=tibble(x=1.5, y=850), fill = NA, label.color = NA, label="*p* < 1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=2, y=1000), fill = NA, label.color = NA, label="*p* = 0.99" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  geom_richtext(data=tibble(x=3.5, y=850), fill = NA, label.color = NA, label="*p* < 1e-7" ,aes(x=x, y=y), inherit.aes=FALSE, size=3)
ggsave(glue("graphs/dbdb_asv_chao.pdf"), height = 6, width = 8, unit = "cm")