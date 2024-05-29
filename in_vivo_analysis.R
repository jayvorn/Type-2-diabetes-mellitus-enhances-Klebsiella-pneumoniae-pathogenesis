library(readxl)
library(tidyverse)
library(glue)
library(dplyr)
library(ggtext)
library(ggplot2)
library(RColorBrewer)
library(DescTools)

# set environment ##################################
setwd("~/Desktop/scratch/dbdb_work")
#===================================================================
# import data ##################################
#Oral gavage intact microbiome: M294-M303 (10), M324-333 (10), M604-613 (9 - no 611), M629-638 (9 - no 638)
M294_303<-read_xlsx("~/Desktop/scratch/dbdb_work/raw_data/20230913_M294_303.xlsx") %>%
  select(Condition, Sample, CFU_per_g_total, LOD, Tissue, Barrier, Genotype) %>%
  filter(Sample != "Inoculum") %>%
  filter(Tissue == "cecum")%>%
  mutate(Abx = "No")
M324_333<-read_xlsx("~/Desktop/scratch/dbdb_work/raw_data/20230920_M324_333.xlsx") %>%
  select(Condition, Sample, CFU_per_g_total, LOD, Tissue, Barrier, Genotype) %>%
  filter(Sample != "Inoculum")%>%
  filter(Tissue == "cecum")%>%
  mutate(Abx = "No")
M604_613<-read_xlsx("~/Desktop/scratch/dbdb_work/raw_data/20240109_M604_613.xlsx") %>%
  select(Condition, Sample, CFU_per_g_total, LOD, Tissue, Barrier, Genotype) %>%
  filter(Sample != "Inoculum") %>%
  filter(Tissue == "cecum")%>%
  filter(Sample != "cecum_M611")%>%
  mutate(Abx = "No")
M629_638<-read_xlsx("~/Desktop/scratch/dbdb_work/raw_data/20240118_M629_638.xlsx") %>%
  select(Condition, Sample, CFU_per_g_total, LOD, Tissue, Barrier, Genotype) %>%
  filter(Sample != "Inoculum") %>%
  filter(Tissue == "cecum")%>%
  filter(Sample != "cecum_M638")%>%
  mutate(Abx = "No")

#Oral gavage disrupted microbiome: M544-553 (10), M689-696 (4 - only 689, 690, 693, 969) 
M544_553<-read_xlsx("~/Desktop/scratch/dbdb_work/raw_data/20231213_M544_553.xlsx") %>%
  select(Condition, Sample, CFU_per_g_total, LOD, Tissue, Barrier, Genotype) %>%
  filter(Sample != "Inoculum")%>%
  filter(Tissue == "cecum")%>%
  mutate(Abx = "Yes")
M689_696<-read_xlsx("~/Desktop/scratch/dbdb_work/raw_data/20240214_M689_698.xlsx") %>%
  select(Condition, Sample, CFU_per_g_total, LOD, Exclude, Tissue, Barrier, Genotype) %>%
  filter(Sample != "Inoculum")%>%
  filter(Tissue == "cecum")%>%
  filter(Exclude == "No") %>%
  select(-Exclude)%>%
  mutate(Abx = "Yes")

#Pneumonia KPPR1 only high dose: M429-448 (20), M699-713 (14 - no 709), M744_757 (14 - includes controls)
M429_438<-read_xlsx("~/Desktop/scratch/dbdb_work/raw_data/20231026_M429_438.xlsx") %>%
  select(Condition, Sample, CFU_per_g_total, CFU_per_mL,LOD, Tissue, Barrier, Genotype) %>%
  filter(Sample != "Inoculum") %>%
  filter(LOD == "Yes" | LOD == "No") %>%
  mutate(dose = "high")
M439_448<-read_xlsx("~/Desktop/scratch/dbdb_work/raw_data/20231102_M439_448.xlsx") %>%
  select(Condition, Sample, CFU_per_g_total, CFU_per_mL,LOD, Tissue, Barrier, Genotype) %>%
  filter(Sample != "Inoculum")  %>%
  filter(LOD == "Yes" | LOD == "No") %>%
  mutate(dose = "high")
M699_713<-read_xlsx("~/Desktop/scratch/dbdb_work/raw_data/20240221_M699_713.xlsx") %>%
  select(Condition, Sample, CFU_per_g_total, CFU_per_mL,LOD, Tissue, Barrier, Genotype) %>%
  filter(Sample != "Inoculum")  %>%
  filter(LOD == "Yes" | LOD == "No") %>%
  mutate(dose = "high")
M744_757<-read_xlsx("~/Desktop/scratch/dbdb_work/raw_data/20240315_M744_757.xlsx") %>%
  select(Condition, Sample, CFU_per_g_total, CFU_per_mL,LOD, Tissue, Barrier, Genotype) %>%
  filter(Sample != "Inoculum")  %>%
  filter(LOD == "Yes" | LOD == "No") %>%
  mutate(dose = "high")

#Pneumonia KPPR1 only low dose: M796_825
M769_825<-read_xlsx("~/Desktop/scratch/dbdb_work/raw_data/20240404_M796_825.xlsx") %>%
  select(Condition, Sample, CFU_per_g_total, CFU_per_mL,LOD, Tissue, Barrier, Genotype) %>%
  filter(Sample != "Inoculum") %>%
  filter(LOD == "Yes" | LOD == "No") %>%
  mutate(dose = "low")

#Pneumonia gltA: M594-603 (10), M614-628 (14 - no 618)
M594_603<-read_xlsx("~/Desktop/scratch/dbdb_work/raw_data/20240109_M594_603.xlsx") %>%
  select(Condition, Sample, Kp_strain_ref,CFU_per_mL_strain, CFU_per_g_strain, LOD, Tissue, Barrier, Genotype, Log_CI) %>%
  filter(Sample != "Inoculum")  %>%
  filter(LOD == "Yes" | LOD == "No")
M614_628<-read_xlsx("~/Desktop/scratch/dbdb_work/raw_data/20240117_M614_628.xlsx") %>%
  select(Condition, Sample, Kp_strain_ref,CFU_per_mL_strain, CFU_per_g_strain, LOD, Tissue, Barrier, Genotype, Log_CI) %>%
  filter(Sample != "Inoculum")  %>%
  filter(LOD == "Yes" | LOD == "No")
#===================================================================
# combine and export data ##################################
#Oral gavage
og_data<-rbind(M294_303, M324_333, M604_613, M629_638, M544_553, M689_696)%>%
  filter(LOD == "No") %>%
  mutate(log_CFU_per_g_total = log(CFU_per_g_total, 10)) %>%
  mutate(group = paste(Genotype, "_", Abx, sep = ""))%>%
  group_by(group)%>%
  mutate(mean_log_CFU_per_g_total = mean(log_CFU_per_g_total)) %>%
  ungroup()
write.csv(og_data, "~/Desktop/scratch/dbdb_work/processed_data/og_data.csv")
og_metadata<-rbind(M294_303, M324_333, M604_613, M629_638, M544_553, M689_696)%>%
  mutate(group = paste(Genotype, "_", Abx, sep = ""))%>%
  mutate(log_CFU_per_g_total = log(CFU_per_g_total, 10)) %>%
  rename(genotype = Genotype) %>%
  rename(barrier = Barrier) %>%
  mutate(sample = paste("M", str_sub(Sample, 7,10), sep = "")) %>%
  select(-Sample) %>%
  select(sample, genotype, barrier, group, log_CFU_per_g_total)
write.csv(og_metadata, "~/Desktop/scratch/dbdb_work/processed_data/og_metadata.csv")
#Lung mono
lung_mono_data<-rbind(M429_438, M439_448, M699_713, M744_757, M769_825)%>%
  filter(LOD == "No")
write.csv(lung_mono_data, "~/Desktop/scratch/dbdb_work/processed_data/lung_mono_data.csv")

lung_mono_data_high<-rbind(M429_438, M439_448, M699_713, M744_757)%>%
  filter(LOD == "No") %>%
  mutate(log_CFU_per_g_total = log(CFU_per_g_total, 10))%>%
  mutate(log_CFU_per_mL = log(CFU_per_mL, 10))%>%
  group_by(Genotype, Tissue)%>%
  mutate(mean_log_CFU_per_g_total = mean(log_CFU_per_g_total)) %>%
  mutate(mean_log_CFU_per_mL = mean(log_CFU_per_mL))%>%
  ungroup()
lung_mono_data_low<-rbind(M769_825)%>%
  filter(LOD == "No") %>%
  mutate(log_CFU_per_g_total = log(CFU_per_g_total, 10))%>%
  mutate(log_CFU_per_mL = log(CFU_per_mL, 10))%>%
  group_by(Genotype, Tissue)%>%
  mutate(mean_log_CFU_per_g_total = mean(log_CFU_per_g_total)) %>%
  mutate(mean_log_CFU_per_mL = mean(log_CFU_per_mL))%>%
  ungroup()
#Lung competative
lung_glta_data_cfu<-rbind(M594_603, M614_628)%>%
  filter(LOD == "No") %>%
  mutate(log_CFU_per_g_strain = log(CFU_per_g_strain, 10))%>%
  mutate(log_CFU_per_mL_strain = log(CFU_per_mL_strain, 10)) %>%
  group_by(Genotype, Tissue, Kp_strain_ref)%>%
  mutate(mean_log_CFU_per_g_strain = mean(log_CFU_per_g_strain))%>%
  mutate(mean_log_CFU_per_mL_strain = mean(log_CFU_per_mL_strain))%>%
  ungroup()
lung_glta_data<-rbind(M594_603, M614_628)%>%
  filter(LOD == "No") %>%
  filter(!is.na(Log_CI)) %>%
  group_by(Genotype, Tissue)%>%
  mutate(mean_Log_CI = mean(Log_CI))%>%
  ungroup()
write.csv(lung_glta_data_cfu, "~/Desktop/scratch/dbdb_work/processed_data/lung_glta_data.csv")
#===================================================================
# define aesthetics ##################################
strain_2_levels<-c("JV1", "JV432")
strain_2_labels<-c("KPPR1", bquote("KPPR1"*Delta*italic(gltA)))
strain_2_colors<-c("black","#50C878")

genotype_3_levels<-c("c57", "het","dbdb")
genotype_3_labels<-c("+/+", "db/+","db/db")
genotype_3_shapes<-c(22,21,24)
genotype_3_colors<-c("black","#5F4B8B","#AD5E99")

genotype_4_levels<-c("het_No","dbdb_No","het_Yes","dbdb_Yes")
genotype_4_labels<-c("db/+<br>No abx","db/db<br>No abx","db/+<br>Abx","db/db<br>Abx")
genotype_4_shapes<-c(21,24, 21,24)
genotype_4_colors<-c("black","black","#5F4B8B","#AD5E99")
genotype_4_fills<-c("#5F4B8B","#AD5E99","grey","lightgrey")
#===================================================================
# oral gavage ##################################
# run stats
og_aov <- TukeyHSD(aov(log_CFU_per_g_total ~ group, data = og_data))
# plot data
ggplot(og_data, aes(x=group, y=log_CFU_per_g_total, color = group, shape = group, fill = group)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=group, y=mean_log_CFU_per_g_total),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  labs(x = NULL, y = bquote(Log~CFU*"/gram cecum")) +
  scale_y_continuous(limits = c(2,13),
                     breaks = c(3,6,9,12)) +
  scale_x_discrete(limits = genotype_4_levels, 
                 labels = genotype_4_labels) +
  scale_color_manual(breaks = genotype_4_levels, 
                     labels = genotype_4_labels,
                     values = genotype_4_colors) +
  scale_fill_manual(breaks = genotype_4_levels, 
                     labels = genotype_4_labels,
                     values = genotype_4_fills) +
  scale_shape_manual(breaks = genotype_4_levels, 
                     labels = genotype_4_labels,
                     values = genotype_4_shapes) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10)) +
  annotate("segment", x = c(1,3), xend = c(2,4), 
         y = c(10, 12.2), yend = c(10,12.2), 
         color = "black") +
  geom_richtext(data=tibble(x=1.5, y=10.4), fill = NA, label.color = NA, label="*p* = 3.6e-3" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=3.5, y=12.6), fill = NA, label.color = NA, label="*p* = 0.99" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=1, y=6.5), fill = NA, label.color = NA, label="*p* = 1.3e-5" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=9.5), fill = NA, label.color = NA, label="*p* = 8.2e-3" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5)
ggsave("graphs/oral_gavage_cfu.pdf", height = 8, width = 6, unit = "cm")
#===================================================================
# lung_mono high dose ##################################
# run stats
lung_mono_data_lung<-lung_mono_data_high %>%
  filter(Tissue == "lung")
lung_mono_aov <- TukeyHSD(aov(log_CFU_per_g_total ~ Genotype, data = lung_mono_data_lung))
liver_mono_data_lung<-lung_mono_data_high %>%
  filter(Tissue == "liver")
liver_mono_aov <- TukeyHSD(aov(log_CFU_per_g_total ~ Genotype, data = liver_mono_data_lung))
blood_mono_data_lung<-lung_mono_data_high %>%
  filter(Tissue == "blood")
blood_mono_aov <- TukeyHSD(aov(log_CFU_per_mL ~ Genotype, data = blood_mono_data_lung))
# plot lung data
ggplot(lung_mono_data_lung, aes(x=Genotype, y=log_CFU_per_g_total, fill = Genotype, shape = Genotype)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=Genotype, y=mean_log_CFU_per_g_total),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  labs(x = NULL, y = bquote(Log~CFU*"/gram lung")) +
  scale_y_continuous(limits = c(5,12),
                     breaks = c(6,8,10,12)) +
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
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10)) +
  annotate("segment", x = c(1,1,2), xend = c(2,3,3), 
           y = c(10,11,7), yend = c(10,11,7), 
           color = "black") +
  geom_richtext(data=tibble(x=1.5, y=10.4), fill = NA, label.color = NA, label="*p* = 0.37" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=11.4), fill = NA, label.color = NA, label="*p* = 0.03" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2.5, y=6.6), fill = NA, label.color = NA, label="*p* = 0.31" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5)
ggsave("graphs/lung_mono_high_cfu.pdf", height = 8, width = 5, unit = "cm")
# plot liver data
ggplot(liver_mono_data_lung, aes(x=Genotype, y=log_CFU_per_g_total, fill = Genotype, shape = Genotype)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=Genotype, y=mean_log_CFU_per_g_total),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  labs(x = NULL, y = bquote(Log~CFU*"/gram liver")) +
  scale_y_continuous(limits = c(2,10),
                     breaks = c(2,4,6,8,10)) +
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
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10)) +
  annotate("segment", x = c(1,1,2), xend = c(2,3,3), 
           y = c(8,9.2,8.2), yend = c(8,9.2,8.2), 
           color = "black") +
  geom_richtext(data=tibble(x=1.5, y=8.4), fill = NA, label.color = NA, label="*p* = 0.67" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=9.6), fill = NA, label.color = NA, label="*p* = 0.41" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2.5, y=8.6), fill = NA, label.color = NA, label="*p* = 0.06" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5)
ggsave("graphs/liver_mono_high_cfu.pdf", height = 8, width = 5, unit = "cm")
# plot blood data
ggplot(blood_mono_data_lung, aes(x=Genotype, y=log_CFU_per_mL, fill = Genotype, shape = Genotype)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=Genotype, y=mean_log_CFU_per_mL),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  labs(x = NULL, y = bquote(Log~CFU*"/mL blood")) +
  scale_y_continuous(limits = c(2,8),
                     breaks = c(2,4,6,8)) +
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
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10)) +
  annotate("segment", x = c(1,1,2), xend = c(2,3,3), 
           y = c(6,7.2,6.2), yend = c(6,7.2,6.2), 
           color = "black") +
  geom_richtext(data=tibble(x=1.5, y=6.4), fill = NA, label.color = NA, label="*p* = 0.89" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=7.6), fill = NA, label.color = NA, label="*p* = 0.27" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2.5, y=6.6), fill = NA, label.color = NA, label="*p* = 0.36" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5)
ggsave("graphs/blood_mono_high_cfu.pdf", height = 8, width = 5, unit = "cm")
#===================================================================
# lung_mono low dose ##################################
# run stats
lung_mono_data_lung<-lung_mono_data_low %>%
  filter(Tissue == "lung")
lung_mono_aov <- TukeyHSD(aov(log_CFU_per_g_total ~ Genotype, data = lung_mono_data_lung))
liver_mono_data_lung<-lung_mono_data_low %>%
  filter(Tissue == "liver")
liver_mono_aov <- TukeyHSD(aov(log_CFU_per_g_total ~ Genotype, data = liver_mono_data_lung))
blood_mono_data_lung<-lung_mono_data_low %>%
  filter(Tissue == "blood")
blood_mono_aov <- TukeyHSD(aov(log_CFU_per_mL ~ Genotype, data = blood_mono_data_lung))
# plot lung data
ggplot(lung_mono_data_lung, aes(x=Genotype, y=log_CFU_per_g_total, fill = Genotype, shape = Genotype)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=Genotype, y=mean_log_CFU_per_g_total),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  labs(x = NULL, y = bquote(Log~CFU*"/gram lung")) +
  scale_y_continuous(limits = c(5,12),
                     breaks = c(6,8,10,12)) +
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
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10)) +
  annotate("segment", x = c(1,1,2), xend = c(2,3,3), 
           y = c(8,9.3,8.5), yend = c(8,9.3,8.5), 
           color = "black") +
  geom_richtext(data=tibble(x=1.5, y=8.4), fill = NA, label.color = NA, label="*p* = 0.36" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=9.7), fill = NA, label.color = NA, label="*p* = 2.4e-4" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2.5, y=8.9), fill = NA, label.color = NA, label="*p* = 6.2e-6" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5)
ggsave("graphs/lung_mono_low_cfu.pdf", height = 8, width = 5, unit = "cm")
# plot liver data
ggplot(liver_mono_data_lung, aes(x=Genotype, y=log_CFU_per_g_total, fill = Genotype, shape = Genotype)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=Genotype, y=mean_log_CFU_per_g_total),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  labs(x = NULL, y = bquote(Log~CFU*"/gram liver")) +
  scale_y_continuous(limits = c(2,10),
                     breaks = c(2,4,6,8,10)) +
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
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10)) +
  annotate("segment", x = c(1,1,2), xend = c(2,3,3), 
           y = c(5,6,5.2), yend = c(5,6,5.2), 
           color = "black") +
  geom_richtext(data=tibble(x=1.5, y=5.4), fill = NA, label.color = NA, label="*p* = 0.65" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=6.4), fill = NA, label.color = NA, label="*p* = 0.98" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2.5, y=5.6), fill = NA, label.color = NA, label="*p* = 0.45" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5)
ggsave("graphs/liver_mono_low_cfu.pdf", height = 8, width = 5, unit = "cm")
# plot blood data
ggplot(blood_mono_data_lung, aes(x=Genotype, y=log_CFU_per_mL, fill = Genotype, shape = Genotype)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=Genotype, y=mean_log_CFU_per_mL),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  labs(x = NULL, y = bquote(Log~CFU*"/mL blood")) +
  scale_y_continuous(limits = c(2,8),
                     breaks = c(2,4,6,8)) +
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
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10)) +
  annotate("segment", x = c(1,1,2), xend = c(2,3,3), 
           y = c(4.5,5.5,4.7), yend = c(4.5,5.5,4.7), 
           color = "black") +
  geom_richtext(data=tibble(x=1.5, y=4.9), fill = NA, label.color = NA, label="*p* = 0.29" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=5.9), fill = NA, label.color = NA, label="*p* = 0.72" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2.5, y=5.1), fill = NA, label.color = NA, label="*p* = 0.59" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5)
ggsave("graphs/blood_mono_low_cfu.pdf", height = 8, width = 5, unit = "cm")
#===================================================================
# lung_glta ##################################
# run stats
lung_glta_data_lung<-lung_glta_data %>%
  filter(Tissue == "lung")
lung_glta_aov <- TukeyHSD(aov(Log_CI ~ Genotype, data = lung_glta_data_lung))
glta_dbdb_lung<-lung_glta_data_lung %>%
  filter(Genotype == "dbdb")
t.test(x = glta_dbdb_lung$Log_CI, rep(0, 9))
glta_het_lung<-lung_glta_data_lung %>%
  filter(Genotype == "het")
t.test(x = glta_het_lung$Log_CI, rep(0, 10))
glta_c57_lung<-lung_glta_data_lung %>%
  filter(Genotype == "c57")
t.test(x = glta_c57_lung$Log_CI, rep(0, 5))
liver_glta_data_lung<-lung_glta_data %>%
  filter(Tissue == "liver") 
liver_glta_aov <- TukeyHSD(aov(Log_CI ~ Genotype, data = liver_glta_data_lung))
glta_dbdb_liver<-liver_glta_data_lung %>%
  filter(Genotype == "dbdb")
t.test(x = glta_dbdb_liver$Log_CI, rep(0, 9))
glta_het_liver<-liver_glta_data_lung %>%
  filter(Genotype == "het")
t.test(x = glta_het_liver$Log_CI, rep(0, 9))
#glta_c57_liver<-liver_glta_data_lung %>%
#  filter(Genotype == "c57")
#t.test(x = glta_c57_liver$Log_CI, rep(0, 1))
blood_glta_data_lung<-lung_glta_data %>%
  filter(Tissue == "blood")
blood_glta_aov <- TukeyHSD(aov(Log_CI ~ Genotype, data = blood_glta_data_lung))
glta_dbdb_blood<-blood_glta_data_lung %>%
  filter(Genotype == "dbdb")
t.test(x = glta_dbdb_blood$Log_CI, rep(0, 8))
glta_het_blood<-blood_glta_data_lung %>%
  filter(Genotype == "het")
t.test(x = glta_het_blood$Log_CI, rep(0, 6))
glta_c57_blood<-blood_glta_data_lung %>%
  filter(Genotype == "c57")
t.test(x = glta_c57_blood$Log_CI, rep(0, 3))
# plot lung data
ggplot(lung_glta_data_lung, aes(x=Genotype, y=Log_CI, fill = Genotype, shape = Genotype)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=Genotype, y=mean_Log_CI),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = bquote(Log~CI["(mutant/WT)"])) +
  scale_y_continuous(limits = c(-4,2),
                     breaks = c(-4,-2,0,2)) +
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
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10)) +
  annotate("segment", x = c(1,1,2), xend = c(2,3,3), 
           y = c(1,1.5,-2.5), yend = c(1,1.5,-2.5), 
           color = "black") +
  geom_richtext(data=tibble(x=1.5, y=1.2), fill = NA, label.color = NA, label="*p* = 4.7e-5" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=1.7), fill = NA, label.color = NA, label="*p* = 1.4e-5" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2.5, y=-2.7), fill = NA, label.color = NA, label="*p* = 0.73" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=1, y=-1.3), fill = NA, label.color = NA, label="*p* = 5.3e-5" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=0.6), fill = NA, label.color = NA, label="*p* = 4.8e-3" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=3, y=0.3), fill = NA, label.color = NA, label="*p* = 2.6e-4" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) 
ggsave("graphs/lung_glta_ci.pdf", height = 8, width = 5, unit = "cm")
# plot liver data
ggplot(liver_glta_data_lung, aes(x=Genotype, y=Log_CI, fill = Genotype, shape = Genotype)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=Genotype, y=mean_Log_CI),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = bquote(Log~CI["(mutant/WT)"])) +
  scale_y_continuous(limits = c(-4,2),
                     breaks = c(-4,-2,0,2)) +
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
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10)) +
  annotate("segment", x = c(2), xend = c(3), 
           y = c(-2), yend = c(-2), 
           color = "black") +
  geom_richtext(data=tibble(x=2.5, y=-2.3), fill = NA, label.color = NA, label="*p* = 0.64" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=2), fill = NA, label.color = NA, label="*p* = 0.2" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=3, y=1.2), fill = NA, label.color = NA, label="*p* = 0.87" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) 
ggsave("graphs/liver_glta_ci.pdf", height = 8, width = 5, unit = "cm")
# plot blood data
ggplot(blood_glta_data_lung, aes(x=Genotype, y=Log_CI, fill = Genotype, shape = Genotype)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=Genotype, y=mean_Log_CI),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = bquote(Log~CI["(mutant/WT)"])) +
  scale_y_continuous(limits = c(-4,2),
                     breaks = c(-4,-2,0,2)) +
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
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10)) +
  annotate("segment", x = c(1,1,2), xend = c(2,3,3), 
           y = c(1,1.5,-2.7), yend = c(1,1.5,-2.7), 
           color = "black") +
  geom_richtext(data=tibble(x=1.5, y=1.2), fill = NA, label.color = NA, label="*p* = 0.99" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=1.7), fill = NA, label.color = NA, label="*p* = 0.81" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2.5, y=-2.9), fill = NA, label.color = NA, label="*p* = 0.64" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=1, y=0.3), fill = NA, label.color = NA, label="*p* = 1.8e-3" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=0.8), fill = NA, label.color = NA, label="*p* = 0.51" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=3, y=1.0), fill = NA, label.color = NA, label="*p* = 0.47" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5)
ggsave("graphs/blood_glta_ci.pdf", height = 8, width = 5, unit = "cm")

# plot lung cfu data
lung_glta_data_cfu %>%
  filter(Tissue == "lung") %>%
  ggplot(aes(x=Genotype, y=log_CFU_per_g_strain, fill = Kp_strain_ref, shape = Genotype, group = Kp_strain_ref)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             size = 2.5,
             show.legend = TRUE) +
  guides(shape = FALSE, fill=guide_legend(override.aes=list(color=strain_2_colors))) +
  geom_point(aes(x=Genotype, y=mean_log_CFU_per_g_strain),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  labs(x = NULL, y = bquote("Log(CFU/gram lung)")) +
  scale_y_continuous(limits = c(1.9,12),
                     breaks = c(3,6,9,12)) +
  scale_x_discrete(limits = genotype_3_levels, 
                   labels = genotype_3_labels) +
  scale_fill_manual(name = "Strain",
                    breaks = strain_2_levels,
                    labels = strain_2_labels,
                    values = strain_2_colors) +
  scale_shape_manual(breaks = genotype_3_levels, 
                     labels = genotype_3_labels,
                     values = genotype_3_shapes) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "right",
        legend.text = element_text(size = 10))
ggsave("graphs/lung_glta_cfu.pdf", height = 8, width = 10, unit = "cm")
# plot liver cfu data
lung_glta_data_cfu %>%
  filter(Tissue == "liver") %>%
  ggplot(aes(x=Genotype, y=log_CFU_per_g_strain, fill = Kp_strain_ref, shape = Genotype, group = Kp_strain_ref)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             size = 2.5,
             show.legend = TRUE) +
  guides(shape = FALSE, fill=guide_legend(override.aes=list(color=strain_2_colors))) +
  geom_point(aes(x=Genotype, y=mean_log_CFU_per_g_strain),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  labs(x = NULL, y = bquote("Log(CFU/gram liver)")) +
  scale_y_continuous(limits = c(1.9,12),
                     breaks = c(3,6,9,12)) +
  scale_x_discrete(limits = genotype_3_levels, 
                   labels = genotype_3_labels) +
  scale_fill_manual(name = "Strain",
                    breaks = strain_2_levels,
                    labels = strain_2_labels,
                    values = strain_2_colors) +
  scale_shape_manual(breaks = genotype_3_levels, 
                     labels = genotype_3_labels,
                     values = genotype_3_shapes) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "right",
        legend.text = element_text(size = 10))
ggsave("graphs/liver_glta_cfu.pdf", height = 8, width = 10, unit = "cm")
# plot liver cfu data
lung_glta_data_cfu %>%
  filter(Tissue == "blood") %>%
  ggplot(aes(x=Genotype, y=log_CFU_per_mL_strain, fill = Kp_strain_ref, shape = Genotype, group = Kp_strain_ref)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             alpha = 1,
             size = 2.5,
             show.legend = TRUE) +
  guides(shape = FALSE, fill=guide_legend(override.aes=list(color=strain_2_colors))) +
  geom_point(aes(x=Genotype, y=mean_log_CFU_per_mL_strain),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  geom_hline(yintercept = log(199, 10), linetype = "dashed", color = "red") +
  labs(x = NULL, y = bquote("Log(CFU/mL blood)")) +
  scale_y_continuous(limits = c(2,12),
                     breaks = c(3,6,9,12)) +
  scale_x_discrete(limits = genotype_3_levels, 
                   labels = genotype_3_labels) +
  scale_fill_manual(name = "Strain",
                    breaks = strain_2_levels,
                    labels = strain_2_labels,
                    values = strain_2_colors) +
  scale_shape_manual(breaks = genotype_3_levels, 
                     labels = genotype_3_labels,
                     values = genotype_3_shapes) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "right",
        legend.text = element_text(size = 10))
ggsave("graphs/blood_glta_cfu.pdf", height = 8, width = 10, unit = "cm")