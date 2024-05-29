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
library(DescTools)
library(gplots)
library(ggrepel)

# set environment ##################################
setwd("~/Desktop/scratch/dbdb_work")
#===================================================================
# define aesthetics ##################################
genotype_3_levels<-c("c57", "het","dbdb")
genotype_3_labels<-c("+/+", "db/+","db/db")
genotype_3_shapes<-c(22,21,24)
genotype_3_colors<-c("black","#5F4B8B","#AD5E99")
#===================================================================
# import data ##################################
M779_795_aa<-read_xlsx("~/Desktop/scratch/dbdb_work/raw_data/AA_batch_4_29_24_final.xlsx", sheet = "processed_data") %>%
  select(-genotype) %>%
  pivot_longer(!sample, names_to ="aa", values_to = "ug_aa_per_ug_protein")
M779_795_key<-read_xlsx("~/Desktop/scratch/dbdb_work/raw_data/AA_batch_4_29_24_final.xlsx", sheet = "processed_data") %>%
  select(sample, genotype)
M779_795_LOD<-read_xlsx("~/Desktop/scratch/dbdb_work/raw_data/AA_batch_4_29_24_final.xlsx", sheet = "raw_data") %>%
  select(`Amino Acid`, LOD_pass) %>%
  rename(aa = `Amino Acid`)
M779_795<-left_join(M779_795_key, M779_795_aa, by = "sample", relationship = "many-to-many") %>%
  left_join(., M779_795_LOD, by = "aa", relationship = "many-to-many")
M779_795 %>%
  filter(LOD_pass > 6) %>%
  select(sample, aa, ug_aa_per_ug_protein)%>%
  pivot_wider(names_from=aa,values_from=ug_aa_per_ug_protein)%>%
  column_to_rownames("sample") %>%
  write.csv(file = "processed_data/balf_aa.csv")
#===================================================================
# run stats ##################################
aa_pvals <- M779_795%>%
  filter(LOD_pass > 6) %>%
  mutate(ug_aa_per_ug_protein = log(ug_aa_per_ug_protein)) %>%
  nest(data = -aa) %>%
  mutate(experiment_tests = map(.x=data, 
                                ~TukeyHSD(aov(ug_aa_per_ug_protein ~ genotype, data=.x)) %>%
                                  tidy())) %>%
  unnest(experiment_tests) %>%
  select(aa, contrast, adj.p.value)

Total_aa<-M779_795%>%
  filter(aa == "Total_aa")
Total_aa_aov<-TukeyHSD(aov(log(ug_aa_per_ug_protein) ~ genotype, data=Total_aa))

aa_pvals_c57_dbdb <- M779_795%>%
  filter(LOD_pass > 6) %>%
  filter(genotype != "het") %>%
  mutate(ug_aa_per_ug_protein = log(ug_aa_per_ug_protein)) %>%
  nest(data = -aa) %>%
  mutate(experiment_tests = map(.x=data, 
                                ~pairwise.t.test(.x$ug_aa_per_ug_protein, .x$genotype, p.adjust.method = "none") %>%
                                  tidy())) %>%
  unnest(experiment_tests) %>%
  select(aa, group1, group2, p.value)

aa_pvals_c57_het <- M779_795%>%
  filter(LOD_pass > 6) %>%
  filter(genotype != "dbdb") %>%
  mutate(ug_aa_per_ug_protein = log(ug_aa_per_ug_protein)) %>%
  nest(data = -aa) %>%
  mutate(experiment_tests = map(.x=data, 
                                ~pairwise.t.test(.x$ug_aa_per_ug_protein, .x$genotype, p.adjust.method = "none") %>%
                                  tidy())) %>%
  unnest(experiment_tests) %>%
  select(aa, group1, group2, p.value)
#================================================
# Create AA dotplots ######################
M779_795 %>%
  filter(aa == "Total_aa") %>%
  group_by(genotype) %>%
  mutate(mean_ug_aa_per_ug_protein = mean(ug_aa_per_ug_protein)) %>%
  ungroup() %>%
  ggplot(., aes(x=genotype, y=ug_aa_per_ug_protein, fill = genotype, shape = genotype)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=genotype, y=mean_ug_aa_per_ug_protein),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  labs(x = NULL, y = bquote(mu*g~"BALF free AAs/"*mu*g~"BALF protein")) +
  scale_y_continuous(limits = c(0, 0.15),
                     breaks = c(0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15)) +
  scale_x_discrete(limits = genotype_3_levels, 
                   labels = genotype_3_labels) +
  scale_fill_manual(breaks = genotype_3_levels, 
                    labels = genotype_3_labels,
                    values = genotype_3_colors) +
  scale_shape_manual(breaks = genotype_3_levels, 
                     labels = genotype_3_labels,
                     values = genotype_3_shapes) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey")) +
  annotate("segment", x = c(1,1), xend = c(2,3), 
           y = c(0.125, 0.14), yend = c(0.125, 0.14), 
           color = "black") +
  geom_richtext(data=tibble(x=1.5, y=0.13), fill = NA, label.color = NA, label="*p* = 0.24" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=0.145), fill = NA, label.color = NA, label="*p* = 0.11" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5)
ggsave("graphs/BALF_AA.pdf", height = 8, width = 5, unit = "cm")

M779_795 %>%
  filter(aa == "Taurine") %>%
  group_by(genotype) %>%
  mutate(mean_ug_aa_per_ug_protein = mean(ug_aa_per_ug_protein)) %>%
  ungroup() %>%
  ggplot(., aes(x=genotype, y=ug_aa_per_ug_protein, fill = genotype, shape = genotype)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=genotype, y=mean_ug_aa_per_ug_protein),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  labs(x = NULL, y = bquote(mu*g~"BALF Taurine/"*mu*g~"BALF protein")) +
  scale_y_continuous(limits = c(0, 0.055),
                     breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05)) +
  scale_x_discrete(limits = genotype_3_levels, 
                   labels = genotype_3_labels) +
  scale_fill_manual(breaks = genotype_3_levels, 
                    labels = genotype_3_labels,
                    values = genotype_3_colors) +
  scale_shape_manual(breaks = genotype_3_levels, 
                     labels = genotype_3_labels,
                     values = genotype_3_shapes) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey")) +
  annotate("segment", x = c(1,1), xend = c(2,3), 
           y = c(0.04, 0.05), yend = c(0.04, 0.05), 
           color = "black") +
  geom_richtext(data=tibble(x=1.5, y=0.042), fill = NA, label.color = NA, label="*p* = 5.3e-3" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=0.052), fill = NA, label.color = NA, label="*p* = 0.02" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5)
ggsave("graphs/BALF_Taurine.pdf", height = 8, width = 5, unit = "cm")

M779_795 %>%
  filter(aa == "Creatine") %>%
  group_by(genotype) %>%
  mutate(mean_ug_aa_per_ug_protein = mean(ug_aa_per_ug_protein)) %>%
  ungroup() %>%
  ggplot(., aes(x=genotype, y=ug_aa_per_ug_protein, fill = genotype, shape = genotype)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=genotype, y=mean_ug_aa_per_ug_protein),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  labs(x = NULL, y = bquote(mu*g~"BALF Creatine/"*mu*g~"BALF protein")) +
  scale_y_continuous(limits = c(0, 0.023),
                     breaks = c(0, 0.005, 0.01, 0.015, 0.02)) +
  scale_x_discrete(limits = genotype_3_levels, 
                   labels = genotype_3_labels) +
  scale_fill_manual(breaks = genotype_3_levels, 
                    labels = genotype_3_labels,
                    values = genotype_3_colors) +
  scale_shape_manual(breaks = genotype_3_levels, 
                     labels = genotype_3_labels,
                     values = genotype_3_shapes) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey")) +
  annotate("segment", x = c(1,1), xend = c(2,3), 
           y = c(0.015, 0.017), yend = c(0.015, 0.017), 
           color = "black") +
  geom_richtext(data=tibble(x=1.5, y=0.0158), fill = NA, label.color = NA, label="*p* = 0.28" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=0.0178), fill = NA, label.color = NA, label="*p* = 0.04" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5)
ggsave("graphs/BALF_Creatine.pdf", height = 8, width = 5, unit = "cm")

M779_795 %>%
  filter(aa == "Glutamine") %>%
  group_by(genotype) %>%
  mutate(mean_ug_aa_per_ug_protein = mean(ug_aa_per_ug_protein)) %>%
  ungroup() %>%
  ggplot(., aes(x=genotype, y=ug_aa_per_ug_protein, fill = genotype, shape = genotype)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=genotype, y=mean_ug_aa_per_ug_protein),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  labs(x = NULL, y = bquote(mu*g~"BALF Glutamine/"*mu*g~"BALF protein")) +
  scale_y_continuous(limits = c(0, 0.012),
                     breaks = c(0, 0.0025, 0.005, 0.0075, 0.01)) +
  scale_x_discrete(limits = genotype_3_levels, 
                   labels = genotype_3_labels) +
  scale_fill_manual(breaks = genotype_3_levels, 
                    labels = genotype_3_labels,
                    values = genotype_3_colors) +
  scale_shape_manual(breaks = genotype_3_levels, 
                     labels = genotype_3_labels,
                     values = genotype_3_shapes) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey")) +
  annotate("segment", x = c(1,1), xend = c(2,3), 
           y = c(0.01, 0.011), yend = c(0.01, 0.011), 
           color = "black") +
  geom_richtext(data=tibble(x=1.5, y=0.0105), fill = NA, label.color = NA, label="*p* = 0.85" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=0.0115), fill = NA, label.color = NA, label="*p* = 0.15" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5)
ggsave("graphs/BALF_Glutamine.pdf", height = 8, width = 5, unit = "cm")

M779_795 %>%
  filter(aa == "Glutamate") %>%
  group_by(genotype) %>%
  mutate(mean_ug_aa_per_ug_protein = mean(ug_aa_per_ug_protein)) %>%
  ungroup() %>%
  ggplot(., aes(x=genotype, y=ug_aa_per_ug_protein, fill = genotype, shape = genotype)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=genotype, y=mean_ug_aa_per_ug_protein),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  labs(x = NULL, y = bquote(mu*g~"BALF Glutamate/"*mu*g~"BALF protein")) +
  scale_y_continuous(limits = c(0, 0.005),
                     breaks = c(0, 0.001, 0.002, 0.003, 0.004, 0.005)) +
  scale_x_discrete(limits = genotype_3_levels, 
                   labels = genotype_3_labels) +
  scale_fill_manual(breaks = genotype_3_levels, 
                    labels = genotype_3_labels,
                    values = genotype_3_colors) +
  scale_shape_manual(breaks = genotype_3_levels, 
                     labels = genotype_3_labels,
                     values = genotype_3_shapes) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey")) +
  annotate("segment", x = c(1,1), xend = c(2,3), 
           y = c(0.0035, 0.004), yend = c(0.0035, 0.004), 
           color = "black") +
  geom_richtext(data=tibble(x=1.5, y=0.0037), fill = NA, label.color = NA, label="*p* = 0.11" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=0.0042), fill = NA, label.color = NA, label="*p* = 0.12" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5)
ggsave("graphs/BALF_Glutamate.pdf", height = 8, width = 5, unit = "cm")

M779_795 %>%
  filter(aa == "Lysine") %>%
  group_by(genotype) %>%
  mutate(mean_ug_aa_per_ug_protein = mean(ug_aa_per_ug_protein)) %>%
  ungroup() %>%
  ggplot(., aes(x=genotype, y=ug_aa_per_ug_protein, fill = genotype, shape = genotype)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=genotype, y=mean_ug_aa_per_ug_protein),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  labs(x = NULL, y = bquote(mu*g~"BALF Lysine/"*mu*g~"BALF protein")) +
  scale_y_continuous(limits = c(0, 0.005),
                     breaks = c(0, 0.001, 0.002, 0.003, 0.004, 0.005)) +
  scale_x_discrete(limits = genotype_3_levels, 
                   labels = genotype_3_labels) +
  scale_fill_manual(breaks = genotype_3_levels, 
                    labels = genotype_3_labels,
                    values = genotype_3_colors) +
  scale_shape_manual(breaks = genotype_3_levels, 
                     labels = genotype_3_labels,
                     values = genotype_3_shapes) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey")) +
  annotate("segment", x = c(1,1), xend = c(2,3), 
           y = c(0.004, 0.0045), yend = c(0.004, 0.0045), 
           color = "black") +
  geom_richtext(data=tibble(x=1.5, y=0.0042), fill = NA, label.color = NA, label="*p* = 0.46" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=0.0047), fill = NA, label.color = NA, label="*p* = 0.92" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5)
ggsave("graphs/BALF_Lysine.pdf", height = 8, width = 5, unit = "cm")
#================================================
# Create AA PCA ######################
# evaluate sample dissimilarity
all_data_pcoa=M779_795 %>%
  filter(LOD_pass > 6) %>%
  filter(aa != "Total_aa" & aa != "gfaa") %>%
  select(sample, aa, ug_aa_per_ug_protein)%>%
  pivot_wider(names_from=aa,values_from=ug_aa_per_ug_protein)%>%
  column_to_rownames("sample")
# generate distance matrix
rda<-rda(all_data_pcoa)

c57_het<-M779_795 %>%
  filter(LOD_pass > 6) %>%
  filter(aa != "Total_aa" & aa != "gfaa") %>%
  filter(genotype != "dbdb") %>%
  select(sample, aa, ug_aa_per_ug_protein) %>%
  pivot_wider(names_from=aa,values_from=ug_aa_per_ug_protein)%>%
  column_to_rownames("sample")
c57_dbdb<-M779_795 %>%
  filter(LOD_pass > 6) %>%
  filter(aa != "Total_aa" & aa != "gfaa") %>%
  filter(genotype != "het") %>%
  select(sample, aa, ug_aa_per_ug_protein) %>%
  pivot_wider(names_from=aa,values_from=ug_aa_per_ug_protein)%>%
  column_to_rownames("sample")
het_dbdb<-M779_795 %>%
  filter(LOD_pass > 6) %>%
  filter(aa != "Total_aa" & aa != "gfaa") %>%
  filter(genotype != "c57") %>%
  select(sample, aa, ug_aa_per_ug_protein) %>%
  pivot_wider(names_from=aa,values_from=ug_aa_per_ug_protein)%>%
  column_to_rownames("sample")

c57_het_metadata<-M779_795_key %>%
  filter(genotype != "dbdb")
c57_dbdb_metadata<-M779_795_key %>%
  filter(genotype != "het")
het_dbdb_metadata<-M779_795_key %>%
  filter(genotype != "c57")

c57_het_dist=vegdist(c57_het, binary=FALSE, method = "euclidean")
c57_dbdb_dist=vegdist(c57_dbdb, binary=FALSE, method = "euclidean")
het_dbdb_dist=vegdist(het_dbdb, binary=FALSE, method = "euclidean")

c57_het_permanova<-adonis2(c57_het_dist ~ c57_het_metadata$genotype, permutations = 1000)
c57_dbdb_permanova<-adonis2(c57_dbdb_dist ~ c57_dbdb_metadata$genotype, permutations = 1000)
het_dbdb_permanova<-adonis2(het_dbdb_dist ~ het_dbdb_metadata$genotype, permutations = 1000)

sample_rda_positions<- as_tibble(rda$CA$u[,1:2], rownames = "sample")
aa_rda_positions<- as_tibble(rda$CA$v[,1:2], rownames = "aa") %>%
  filter(abs(PC1) > 0.09 | abs(PC2) > 0.09) 

sample_rda_positions %>%
  left_join(., M779_795_key, by = "sample") %>%
  ggplot(aes(x=PC1, y=PC2, fill=genotype, shape = genotype)) +
  geom_point(alpha = 1, size = 2.5) +
  stat_ellipse(level=0.7,
               geom="polygon",
               alpha=0.1,
               show.legend=FALSE) +
  geom_segment(data = aa_rda_positions, 
               aes(x = 0, xend = PC1, y = 0, yend = PC2), 
               color = "red",
               arrow = arrow(type = "closed", length = unit(0.1, "cm")),
               inherit.aes = FALSE) +
  geom_text_repel(data = aa_rda_positions, 
                   aes(x = PC1, y = PC2, label = aa), 
                   color = "black",
                  size = 3,
                  box.padding = 2, 
                  inherit.aes = FALSE) +
  labs(x = "Axis 1 (92.9%)", y = "Axis 2 (4.4%)") +
  scale_fill_manual(name = "Genotype",
                     breaks = genotype_3_levels,
                     values = genotype_3_colors,
                     labels = genotype_3_labels) +
  scale_shape_manual(name = "Genotype",
                     breaks = genotype_3_levels,
                     values = genotype_3_shapes,
                     labels = genotype_3_labels)+
  scale_x_continuous(limits = c(-1, 0.7),
                     breaks = c(-1, -0.5, 0, 0.5),
                     labels = c(-1, -0.5, 0, 0.5)) +
  scale_y_continuous(limits = c(-0.8, 1),
                     breaks = c(-0.5, 0, 0.5, 1),
                     labels = c(-0.5, 0, 0.5, 1)) +
  theme_classic() +
  theme(axis.title = element_markdown(size = 14),
        axis.text.x = element_markdown(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10)) +
  guides(color = guide_legend(nrow = 2))
ggsave("graphs/aa_PCoA.pdf", height = 10, width = 10, unit = "cm")
