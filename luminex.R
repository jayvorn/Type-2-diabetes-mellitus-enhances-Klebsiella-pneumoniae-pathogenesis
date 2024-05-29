library(readxl)
library(tidyverse)
library(glue)
library(dplyr)
library(ggtext)
library(ggplot2)
library(RColorBrewer)
library(DescTools)
library(Hmisc)

# set environment ##################################
setwd("~/Desktop/scratch/dbdb_work")
#===================================================================
# define aesthetics ##################################
breaks<-c("Sham_c57", "Sham_het", "Sham_dbdb","Exp_c57", "Exp_het", "Exp_dbdb")
labels<-c("+/+", "db/+", "db/db","+/+", "db/+", "db/db")
shapes<-c(22,21,24,22,21,24)
fills<-c("black","#5F4B8B","#AD5E99","black","#5F4B8B","#AD5E99")
color_levels<-c("Yes", "No")
colors<-c("red","black")
#===================================================================
# import data #################################
#Pneumonia KPPR1 only: M429-448 (20), M699-713 (14 - no 709), M744_757 (14 - includes controls)
#Lung mono
lung_mono_data<-read.csv(lung_mono_data, "~/Desktop/scratch/dbdb_work/processed_data/lung_mono_data.csv")
# luminex data
luminex_data<-read_xlsx("~/Desktop/scratch/dbdb_work/raw_data/20240315_lung_luminex_final.xlsx") %>%
  rename(Sample = 'Sample ID') 
comb_data<-left_join(lung_mono_data, luminex_data, by = "Sample")%>%
  mutate(Group = case_when(grepl("sham", Condition) ~ "Sham",
                           TRUE ~ "Exp"))%>%
  mutate(Group_genotype = paste(Group, Genotype, sep = "_"))%>%
  mutate(log_CFU_per_g_total = log(CFU_per_g_total)) %>%
  mutate(IFN_pg_per_g = ((IFN_pg_ml/3)/Sample_weight_mg)*1000)%>%
  mutate(log_IFN_pg_per_g = log(IFN_pg_per_g))%>%
  mutate(IFN_LOD = case_when(IFN_pg_ml == 1 ~ "Yes",
                             TRUE ~ "No"))%>%
  group_by(Group_genotype)%>%
  mutate(mean_IFN_pg_per_g = mean(IFN_pg_per_g))%>%
  ungroup()%>%
  mutate(IL1_pg_per_g = ((IL1_pg_ml/3)/Sample_weight_mg)*1000)%>%
  mutate(log_IL1_pg_per_g = log(IL1_pg_per_g))%>%
  mutate(IL1_LOD = case_when(IL1_pg_ml == 1 ~ "Yes",
                             TRUE ~ "No"))%>%
  group_by(Group_genotype)%>%
  mutate(mean_IL1_pg_per_g = mean(IL1_pg_per_g))%>%
  ungroup()%>%
  mutate(IL6_pg_per_g = ((IL6_pg_ml/3)/Sample_weight_mg)*1000)%>%
  mutate(log_IL6_pg_per_g = log(IL6_pg_per_g))%>%
  mutate(IL6_LOD = case_when(IL6_pg_ml == 1 ~ "Yes",
                             IL6_pg_ml == 30507.31 ~ "Yes",
                             TRUE ~ "No"))%>%
  group_by(Group_genotype)%>%
  mutate(mean_IL6_pg_per_g = mean(IL6_pg_per_g))%>%
  ungroup()%>%
  mutate(IL10_pg_per_g = ((IL10_pg_ml/3)/Sample_weight_mg)*1000)%>%
  mutate(log_IL10_pg_per_g = log(IL10_pg_per_g))%>%
  mutate(IL10_LOD = case_when(IL10_pg_ml == 1 ~ "Yes",
                             TRUE ~ "No"))%>%
  group_by(Group_genotype)%>%
  mutate(mean_IL10_pg_per_g = mean(IL10_pg_per_g))%>%
  ungroup()%>%
  mutate(IL17_pg_per_g = ((IL17_pg_ml/3)/Sample_weight_mg)*1000)%>%
  mutate(log_IL17_pg_per_g = log(IL17_pg_per_g))%>%
  mutate(IL17_LOD = case_when(IL17_pg_ml == 1 ~ "Yes",
                             TRUE ~ "No"))%>%
  group_by(Group_genotype)%>%
  mutate(mean_IL17_pg_per_g = mean(IL17_pg_per_g))%>%
  ungroup()%>%
  mutate(CCL2_pg_per_g = ((CCL2_pg_ml/3)/Sample_weight_mg)*1000)%>%
  mutate(log_CCL2_pg_per_g = log(CCL2_pg_per_g))%>%
  mutate(CCL2_LOD = case_when(CCL2_pg_ml == 1 ~ "Yes",
                             TRUE ~ "No"))%>%
  group_by(Group_genotype)%>%
  mutate(mean_CCL2_pg_per_g = mean(CCL2_pg_per_g))%>%
  ungroup()%>%
  mutate(TNF_pg_per_g = ((TNF_pg_ml/3)/Sample_weight_mg)*1000)%>%
  mutate(log_TNF_pg_per_g = log(TNF_pg_per_g))%>%
  mutate(TNF_LOD = case_when(TNF_pg_ml == 1 ~ "Yes",
                             TRUE ~ "No"))%>%
  group_by(Group_genotype)%>%
  mutate(mean_TNF_pg_per_g = mean(TNF_pg_per_g))%>%
  ungroup()%>%
  mutate(VEGF_pg_per_g = ((VEGF_pg_ml/3)/Sample_weight_mg)*1000)%>%
  mutate(log_VEGF_pg_per_g = log(VEGF_pg_per_g))%>%
  mutate(VEGF_LOD = case_when(VEGF_pg_ml == 1 ~ "Yes",
                             TRUE ~ "No")) %>%
  group_by(Group_genotype)%>%
  mutate(mean_VEGF_pg_per_g = mean(VEGF_pg_per_g))%>%
  ungroup()%>%
  select(-Sample, Tissue, Sample_weight_mg)
#===================================================================
# run stats #################################
IFN_exp_stat<-comb_data %>%
  filter(Group == "Exp") %>%
  select(Genotype, log_CFU_per_g_total, log_IFN_pg_per_g, IFN_LOD) %>%
  filter(IFN_LOD == "No")
IFN_exp_aov<-TukeyHSD(aov(log_IFN_pg_per_g ~ Genotype, data = IFN_exp_stat))
IFN_sham_stat<-comb_data %>%
  filter(Group == "Sham") %>%
  select(Genotype, log_IFN_pg_per_g, IFN_LOD) %>%
  filter(IFN_LOD == "No")
IFN_sham_aov<-TukeyHSD(aov(log_IFN_pg_per_g ~ Genotype, data = IFN_sham_stat))

IL1_exp_stat<-comb_data %>%
  filter(Group == "Exp") %>%
  select(Genotype, log_CFU_per_g_total, log_IL1_pg_per_g, IL1_LOD) %>%
  filter(IL1_LOD == "No")
IL1_exp_aov<-TukeyHSD(aov(log_IL1_pg_per_g ~ Genotype, data = IL1_exp_stat))
IL1_sham_stat<-comb_data %>%
  filter(Group == "Sham") %>%
  select(Genotype, log_IL1_pg_per_g, IL1_LOD) %>%
  filter(IL1_LOD == "No")
IL1_sham_aov<-TukeyHSD(aov(log_IL1_pg_per_g ~ Genotype, data = IL1_sham_stat))

IL6_exp_stat<-comb_data %>%
  filter(Group == "Exp") %>%
  select(Genotype, log_CFU_per_g_total, log_IL6_pg_per_g, IL6_LOD) %>%
  filter(IL6_LOD == "No")
IL6_exp_aov<-TukeyHSD(aov(log_IL6_pg_per_g ~ Genotype, data = IL6_exp_stat))
IL6_sham_stat<-comb_data %>%
  filter(Group == "Sham") %>%
  select(Genotype, log_IL6_pg_per_g, IL6_LOD) %>%
  filter(IL6_LOD == "No")
IL6_sham_aov<-TukeyHSD(aov(log_IL6_pg_per_g ~ Genotype, data = IL6_sham_stat))

IL10_exp_stat<-comb_data %>%
  filter(Group == "Exp") %>%
  select(Genotype, log_CFU_per_g_total, log_IL10_pg_per_g, IL10_LOD) %>%
  filter(IL10_LOD == "No")
IL10_exp_aov<-TukeyHSD(aov(log_IL10_pg_per_g ~ Genotype, data = IL10_exp_stat))
IL10_sham_stat<-comb_data %>%
  filter(Group == "Sham") %>%
  select(Genotype, log_IL10_pg_per_g, IL10_LOD) %>%
  filter(IL10_LOD == "No")
IL10_sham_aov<-TukeyHSD(aov(log_IL10_pg_per_g ~ Genotype, data = IL10_sham_stat))

IL17_exp_stat<-comb_data %>%
  filter(Group == "Exp") %>%
  select(Genotype, log_CFU_per_g_total, log_IL17_pg_per_g, IL17_LOD) %>%
  filter(IL17_LOD == "No")
IL17_exp_aov<-TukeyHSD(aov(log_IL17_pg_per_g ~ Genotype, data = IL17_exp_stat))
IL17_sham_stat<-comb_data %>%
  filter(Group == "Sham") %>%
  select(Genotype, log_IL17_pg_per_g, IL17_LOD) %>%
  filter(IL17_LOD == "No")
IL17_sham_aov<-TukeyHSD(aov(log_IL17_pg_per_g ~ Genotype, data = IL17_sham_stat))

CCL2_exp_stat<-comb_data %>%
  filter(Group == "Exp") %>%
  select(Genotype, log_CFU_per_g_total, log_CCL2_pg_per_g, CCL2_LOD) %>%
  filter(CCL2_LOD == "No")
CCL2_exp_aov<-TukeyHSD(aov(log_CCL2_pg_per_g ~ Genotype, data = CCL2_exp_stat))
CCL2_sham_stat<-comb_data %>%
  filter(Group == "Sham") %>%
  select(Genotype, log_CCL2_pg_per_g, CCL2_LOD) %>%
  filter(CCL2_LOD == "No")
CCL2_sham_aov<-TukeyHSD(aov(log_CCL2_pg_per_g ~ Genotype, data = CCL2_sham_stat))

TNF_exp_stat<-comb_data %>%
  filter(Group == "Exp") %>%
  select(Genotype,log_CFU_per_g_total, log_TNF_pg_per_g, TNF_LOD) %>%
  filter(TNF_LOD == "No")
TNF_exp_aov<-TukeyHSD(aov(log_TNF_pg_per_g ~ Genotype, data = TNF_exp_stat))
TNF_sham_stat<-comb_data %>%
  filter(Group == "Sham") %>%
  select(Genotype, log_TNF_pg_per_g, TNF_LOD) %>%
  filter(TNF_LOD == "No")
TNF_sham_aov<-TukeyHSD(aov(log_TNF_pg_per_g ~ Genotype, data = TNF_sham_stat))

VEGF_exp_stat<-comb_data %>%
  filter(Group == "Exp") %>%
  select(Genotype, log_CFU_per_g_total, log_VEGF_pg_per_g, VEGF_LOD) %>%
  filter(VEGF_LOD == "No")
VEGF_exp_aov<-TukeyHSD(aov(log_VEGF_pg_per_g ~ Genotype, data = VEGF_exp_stat))
VEGF_sham_stat<-comb_data %>%
  filter(Group == "Sham") %>%
  select(Genotype, log_VEGF_pg_per_g, VEGF_LOD) %>%
  filter(VEGF_LOD == "No")
VEGF_sham_aov<-TukeyHSD(aov(log_VEGF_pg_per_g ~ Genotype, data = VEGF_sham_stat))
#===================================================================
# graph luminex results #################################
#IFN
comb_data %>%
  ggplot(aes(x = Group_genotype, y = IFN_pg_per_g, fill = Group_genotype, shape = Group_genotype)) +
  geom_point(aes(color = IFN_LOD, group = Group_genotype), 
             position = position_jitterdodge(jitter.width = 0.2),
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=Group_genotype, y=mean_IFN_pg_per_g),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  geom_vline(xintercept = 3.5, color = "black", linetype = "dashed") +
  geom_richtext(data=tibble(x=1, y=50), fill = NA, label.color = NA, label="Sham",aes(x=x, y=y), inherit.aes=FALSE, size=3.5) +
  geom_richtext(data=tibble(x=4, y=50), fill = NA, label.color = NA, label="Kp",aes(x=x, y=y), inherit.aes=FALSE, size=3.5) +
  labs(x = NULL, y = bquote(pg~IFN*gamma*"/gram lung")) +
  scale_x_discrete(limits = breaks,
                   breaks = breaks,
                   labels = labels) +
  scale_y_continuous(limits = c(0,50),
                   breaks = c(0,10,20,30,40,50),
                   labels = c(0,10,20,30,40,50))+
  scale_shape_manual(breaks = breaks,
                     values = shapes,
                     labels = labels) +
  scale_fill_manual(breaks = breaks,
                    values = fills,
                    labels = labels) +
  scale_color_manual(breaks = color_levels,
                     values = colors) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10)) +
  annotate("segment", x = c(4,4,5), xend = c(5,6,6), 
           y = c(35,42,37), yend = c(35,42,37), 
           color = "black") +
  geom_richtext(data=tibble(x=4.5, y=37), fill = NA, label.color = NA, label="*p* = 0.85" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=5, y=44), fill = NA, label.color = NA, label="*p* = 0.34" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=5.5, y=39), fill = NA, label.color = NA, label="*p* = 0.57" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5)
ggsave("graphs/lung_IFN.pdf", height = 6, width = 8, unit = "cm")
#IL1
comb_data %>%
  ggplot(aes(x = Group_genotype, y = IL1_pg_per_g, fill = Group_genotype, shape = Group_genotype)) +
  geom_point(aes(color = IL1_LOD, group = Group_genotype), 
             position = position_jitterdodge(jitter.width = 0.2),
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=Group_genotype, y=mean_IL1_pg_per_g),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  geom_vline(xintercept = 3.5, color = "black", linetype = "dashed") +
  geom_richtext(data=tibble(x=1, y=800), fill = NA, label.color = NA, label="Sham",aes(x=x, y=y), inherit.aes=FALSE, size=3.5) +
  geom_richtext(data=tibble(x=4, y=800), fill = NA, label.color = NA, label="Kp",aes(x=x, y=y), inherit.aes=FALSE, size=3.5) +
  labs(x = NULL, y = bquote(pg~IL*"-1"*beta*"/gram lung")) +
  scale_x_discrete(limits = breaks,
                   breaks = breaks,
                   labels = labels) +
  scale_y_continuous(limits = c(0,800),
                     breaks = c(0,200,400,600,800),
                     labels = c(0,200,400,600,800))+
  scale_shape_manual(breaks = breaks,
                     values = shapes,
                     labels = labels) +
  scale_fill_manual(breaks = breaks,
                    values = fills,
                    labels = labels) +
  scale_color_manual(breaks = color_levels,
                     values = colors) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10)) +
  annotate("segment", x = c(4,4,5,1,1,2), xend = c(5,6,6,2,3,3), 
           y = c(400,720,650,50,200,120), yend = c(400,720,650,50,200,120), 
           color = "black") +
  geom_richtext(data=tibble(x=4.5, y=430), fill = NA, label.color = NA, label="*p* = 2.6e-3" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=5, y=750), fill = NA, label.color = NA, label="*p* = 2.7e-6" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=5.5, y=680), fill = NA, label.color = NA, label="*p* = 0.03" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=1.5, y=80), fill = NA, label.color = NA, label="*p* = 0.12" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=230), fill = NA, label.color = NA, label="*p* = 0.37" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2.5, y=150), fill = NA, label.color = NA, label="*p* = 0.65" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) 
ggsave("graphs/lung_IL1.pdf", height = 6, width = 8, unit = "cm")
#IL6
comb_data %>%
  ggplot(aes(x = Group_genotype, y = IL6_pg_per_g, fill = Group_genotype, shape = Group_genotype)) +
  geom_point(aes(color = IL6_LOD, group = Group_genotype), 
             position = position_jitterdodge(jitter.width = 0.2),
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=Group_genotype, y=mean_IL6_pg_per_g),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  geom_vline(xintercept = 3.5, color = "black", linetype = "dashed") +
  geom_richtext(data=tibble(x=1, y=50000), fill = NA, label.color = NA, label="Sham",aes(x=x, y=y), inherit.aes=FALSE, size=3.5) +
  geom_richtext(data=tibble(x=4, y=50000), fill = NA, label.color = NA, label="Kp",aes(x=x, y=y), inherit.aes=FALSE, size=3.5) +
  labs(x = NULL, y = bquote(pg~IL*"-6"*"/gram lung")) +
  scale_x_discrete(limits = breaks,
                   breaks = breaks,
                   labels = labels) +
  scale_y_continuous(limits = c(0,50000),
                     breaks = c(0,10000,20000,30000,40000,50000),
                     labels = c(0,10000,20000,30000,40000,50000))+
  scale_shape_manual(breaks = breaks,
                     values = shapes,
                     labels = labels) +
  scale_fill_manual(breaks = breaks,
                    values = fills,
                    labels = labels) +
  scale_color_manual(breaks = color_levels,
                     values = colors) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10)) +
  annotate("segment", x = c(4,4,5,1,1,2), xend = c(5,6,6,2,3,3), 
           y = c(15000,45000,40000,2000,10000,5000), yend = c(15000,45000,40000,2000,10000,5000), 
           color = "black") +
  geom_richtext(data=tibble(x=4.5, y=17000), fill = NA, label.color = NA, label="*p* = 0.36" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=5, y=47000), fill = NA, label.color = NA, label="*p* = 0.66" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=5.5, y=42000), fill = NA, label.color = NA, label="*p* = 0.82" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=1.5, y=4000), fill = NA, label.color = NA, label="*p* = 0.48" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=12000), fill = NA, label.color = NA, label="*p* = 0.94" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2.5, y=7000), fill = NA, label.color = NA, label="*p* = 0.74" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) 
ggsave("graphs/lung_IL6.pdf", height = 6, width = 8, unit = "cm")
#IL10
comb_data %>%
  ggplot(aes(x = Group_genotype, y = IL10_pg_per_g, fill = Group_genotype, shape = Group_genotype)) +
  geom_point(aes(color = IL10_LOD, group = Group_genotype), 
             position = position_jitterdodge(jitter.width = 0.2),
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=Group_genotype, y=mean_IL10_pg_per_g),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  geom_vline(xintercept = 3.5, color = "black", linetype = "dashed") +
  geom_richtext(data=tibble(x=1, y=350), fill = NA, label.color = NA, label="Sham",aes(x=x, y=y), inherit.aes=FALSE, size=3.5) +
  geom_richtext(data=tibble(x=4, y=350), fill = NA, label.color = NA, label="Kp",aes(x=x, y=y), inherit.aes=FALSE, size=3.5) +
  labs(x = NULL, y = bquote(pg~IL*"-10"*"/gram lung")) +
  scale_x_discrete(limits = breaks,
                   breaks = breaks,
                   labels = labels) +
  scale_y_continuous(limits = c(0,350),
                     breaks = c(0,100,200,300),
                     labels = c(0,100,200,300))+
  scale_shape_manual(breaks = breaks,
                     values = shapes,
                     labels = labels) +
  scale_fill_manual(breaks = breaks,
                    values = fills,
                    labels = labels) +
  scale_color_manual(breaks = color_levels,
                     values = colors) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10)) +
  annotate("segment", x = c(4,4,5,1,1,2), xend = c(5,6,6,2,3,3), 
           y = c(250,310,270,50,110,70), yend = c(250,310,270,50,110,70), 
           color = "black") +
  geom_richtext(data=tibble(x=4.5, y=265), fill = NA, label.color = NA, label="*p* = 0.21" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=5, y=325), fill = NA, label.color = NA, label="*p* = 0.94" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=5.5, y=285), fill = NA, label.color = NA, label="*p* = 0.25" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=1.5, y=65), fill = NA, label.color = NA, label="*p* = 0.48" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=125), fill = NA, label.color = NA, label="*p* = 0.68" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2.5, y=85), fill = NA, label.color = NA, label="*p* = 0.26" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) 
ggsave("graphs/lung_IL10.pdf", height = 6, width = 8, unit = "cm")
#IL17
comb_data %>%
  ggplot(aes(x = Group_genotype, y = IL17_pg_per_g, fill = Group_genotype, shape = Group_genotype)) +
  geom_point(aes(color = IL17_LOD, group = Group_genotype), 
             position = position_jitterdodge(jitter.width = 0.2),
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=Group_genotype, y=mean_IL17_pg_per_g),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  geom_vline(xintercept = 3.5, color = "black", linetype = "dashed") +
  geom_richtext(data=tibble(x=1, y=500), fill = NA, label.color = NA, label="Sham",aes(x=x, y=y), inherit.aes=FALSE, size=3.5) +
  geom_richtext(data=tibble(x=4, y=500), fill = NA, label.color = NA, label="Kp",aes(x=x, y=y), inherit.aes=FALSE, size=3.5) +
  labs(x = NULL, y = bquote("pg IL-17A/gram lung")) +
  scale_x_discrete(limits = breaks,
                   breaks = breaks,
                   labels = labels) +
  scale_y_continuous(limits = c(0,500),
                     breaks = c(0,100,200,300,400,500),
                     labels = c(0,100,200,300,400,500))+
  scale_shape_manual(breaks = breaks,
                     values = shapes,
                     labels = labels) +
  scale_fill_manual(breaks = breaks,
                    values = fills,
                    labels = labels) +
  scale_color_manual(breaks = color_levels,
                     values = colors) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10)) +
  annotate("segment", x = c(4,4,5), xend = c(5,6,6), 
           y = c(300,400,350), yend = c(300,400,350), 
           color = "black") +
  geom_richtext(data=tibble(x=4.5, y=320), fill = NA, label.color = NA, label="*p* = 0.99" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=5, y=420), fill = NA, label.color = NA, label="*p* = 0.98" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=5.5, y=370), fill = NA, label.color = NA, label="*p* = 0.97" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5)
ggsave("graphs/lung_IL17.pdf", height = 6, width = 8, unit = "cm")
#CCL2
comb_data %>%
  ggplot(aes(x = Group_genotype, y = CCL2_pg_per_g, fill = Group_genotype, shape = Group_genotype)) +
  geom_point(aes(color = CCL2_LOD, group = Group_genotype), 
             position = position_jitterdodge(jitter.width = 0.2),
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=Group_genotype, y=mean_CCL2_pg_per_g),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  geom_vline(xintercept = 3.5, color = "black", linetype = "dashed") +
  geom_richtext(data=tibble(x=1, y=7000), fill = NA, label.color = NA, label="Sham",aes(x=x, y=y), inherit.aes=FALSE, size=3.5) +
  geom_richtext(data=tibble(x=4, y=7000), fill = NA, label.color = NA, label="Kp",aes(x=x, y=y), inherit.aes=FALSE, size=3.5) +
  labs(x = NULL, y = bquote("pg CCL2/gram lung")) +
  scale_x_discrete(limits = breaks,
                   breaks = breaks,
                   labels = labels) +
  scale_y_continuous(limits = c(0,7000),
                     breaks = c(0,2000,4000,6000),
                     labels = c(0,2000,4000,6000))+
  scale_shape_manual(breaks = breaks,
                     values = shapes,
                     labels = labels) +
  scale_fill_manual(breaks = breaks,
                    values = fills,
                    labels = labels) +
  scale_color_manual(breaks = color_levels,
                     values = colors) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10)) +
  annotate("segment", x = c(4,4,5), xend = c(5,6,6), 
           y = c(3200,6000,4000), yend = c(3200,6000,4000), 
           color = "black") +
  geom_richtext(data=tibble(x=4.5, y=3500), fill = NA, label.color = NA, label="*p* = 6.4e-4" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=5, y=6300), fill = NA, label.color = NA, label="*p* = 0.052" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=5.5, y=4300), fill = NA, label.color = NA, label="*p* = 0.16" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5)
ggsave("graphs/lung_CCL2.pdf", height = 6, width = 8, unit = "cm")
#TNF
comb_data %>%
  ggplot(aes(x = Group_genotype, y = TNF_pg_per_g, fill = Group_genotype, shape = Group_genotype)) +
  geom_point(aes(color = TNF_LOD, group = Group_genotype), 
             position = position_jitterdodge(jitter.width = 0.2),
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=Group_genotype, y=mean_TNF_pg_per_g),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  geom_vline(xintercept = 3.5, color = "black", linetype = "dashed") +
  geom_richtext(data=tibble(x=1, y=4000), fill = NA, label.color = NA, label="Sham",aes(x=x, y=y), inherit.aes=FALSE, size=3.5) +
  geom_richtext(data=tibble(x=4, y=4000), fill = NA, label.color = NA, label="Kp",aes(x=x, y=y), inherit.aes=FALSE, size=3.5) +
  labs(x = NULL, y = bquote(pg~TNF*alpha*"/gram lung")) +
  scale_x_discrete(limits = breaks,
                   breaks = breaks,
                   labels = labels) +
  scale_y_continuous(limits = c(0,4000),
                     breaks = c(0,1000,2000,3000,4000),
                     labels = c(0,1000,2000,3000,4000))+
  scale_shape_manual(breaks = breaks,
                     values = shapes,
                     labels = labels) +
  scale_fill_manual(breaks = breaks,
                    values = fills,
                    labels = labels) +
  scale_color_manual(breaks = color_levels,
                     values = colors) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10)) +
  annotate("segment", x = c(4,4,5,1,1,2), xend = c(5,6,6,2,3,3), 
           y = c(2500,3300,2900,300,1000,600), yend = c(2500,3300,2900,300,1000,600), 
           color = "black") +
  geom_richtext(data=tibble(x=4.5, y=2650), fill = NA, label.color = NA, label="*p* = 7.2e-4" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=5, y=3450), fill = NA, label.color = NA, label="*p* = 2.3e-3" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=5.5, y=3050), fill = NA, label.color = NA, label="*p* = 0.91" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=1.5, y=450), fill = NA, label.color = NA, label="*p* = 0.28" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=1150), fill = NA, label.color = NA, label="*p* = 0.64" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2.5, y=750), fill = NA, label.color = NA, label="*p* = 0.09" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) 
ggsave("graphs/lung_TNF.pdf", height = 6, width = 8, unit = "cm")
#VEGF
comb_data %>%
  ggplot(aes(x = Group_genotype, y = VEGF_pg_per_g, fill = Group_genotype, shape = Group_genotype)) +
  geom_point(aes(color = VEGF_LOD, group = Group_genotype), 
             position = position_jitterdodge(jitter.width = 0.2),
             size = 2.5,
             show.legend = FALSE) +
  geom_point(aes(x=Group_genotype, y=mean_VEGF_pg_per_g),
             position = position_dodge(width = 0.75),
             color = "black",
             shape = 95,
             size = 8,
             show.legend = FALSE) +
  geom_vline(xintercept = 3.5, color = "black", linetype = "dashed") +
  geom_richtext(data=tibble(x=1, y=600), fill = NA, label.color = NA, label="Sham",aes(x=x, y=y), inherit.aes=FALSE, size=3.5) +
  geom_richtext(data=tibble(x=4, y=600), fill = NA, label.color = NA, label="Kp",aes(x=x, y=y), inherit.aes=FALSE, size=3.5) +
  labs(x = NULL, y = bquote("pg VEGF-A/gram lung")) +
  scale_x_discrete(limits = breaks,
                   breaks = breaks,
                   labels = labels) +
  scale_y_continuous(limits = c(0,600),
                     breaks = c(0,150,300,450,600),
                     labels = c(0,150,300,450,600))+
  scale_shape_manual(breaks = breaks,
                     values = shapes,
                     labels = labels) +
  scale_fill_manual(breaks = breaks,
                    values = fills,
                    labels = labels) +
  scale_color_manual(breaks = color_levels,
                     values = colors) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10)) +
  annotate("segment", x = c(4,4,5,1,1,2), xend = c(5,6,6,2,3,3), 
           y = c(350,550,490,150,375,310), yend = c(350,550,490,150,375,310), 
           color = "black") +
  geom_richtext(data=tibble(x=4.5, y=375), fill = NA, label.color = NA, label="*p* = 6.1e-5" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=5, y=575), fill = NA, label.color = NA, label="*p* < 1e-15" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=5.5, y=515), fill = NA, label.color = NA, label="*p* = 4.3e-5" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=1.5, y=175), fill = NA, label.color = NA, label="*p* = 0.49" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2, y=400), fill = NA, label.color = NA, label="*p* = 5.8e-4" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=2.5, y=335), fill = NA, label.color = NA, label="*p* = 1.4e-3" ,aes(x=x, y=y), inherit.aes=FALSE, size=2.5) 
ggsave("graphs/lung_VEGF.pdf", height = 6, width = 8, unit = "cm")
#===================================================================
# correlate luminex results #################################
#IL1
rcorr(IL1_exp_stat$log_CFU_per_g_total, IL1_exp_stat$log_IL1_pg_per_g)
ggplot(IL1_exp_stat, aes(x = log_IL1_pg_per_g, y = log_CFU_per_g_total, fill = Genotype, shape = Genotype)) +
  geom_point(size = 2.5,show.legend = FALSE) +
  geom_richtext(data=tibble(x=4.25, y=22), fill = NA, label.color = NA, label="r = 0.62",aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=4.25, y=21.5), fill = NA, label.color = NA, label="*p* = 0.00",aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  labs(x = "Log(CFU/g lung)", y = bquote("Log("*IL1*beta*"/gram lung)")) +
  scale_shape_manual(breaks = c("c57", "het", "dbdb"),
                     values = c(22,21,21),
                     labels = c("+/+", "db/+", "db/db")) +
  scale_fill_manual(breaks = c("c57", "het", "dbdb"),
                    values = c("black","#5F4B8B","#AD5E99"),
                    labels = c("+/+", "db/+", "db/db")) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10))
ggsave("graphs/lung_IL1_corr.pdf", height = 6, width = 6, unit = "cm")
#CCL2
rcorr(CCL2_exp_stat$log_CFU_per_g_total, CCL2_exp_stat$log_CCL2_pg_per_g)
ggplot(CCL2_exp_stat, aes(x = log_CCL2_pg_per_g, y = log_CFU_per_g_total, fill = Genotype, shape = Genotype)) +
  geom_point(size = 2.5,show.legend = FALSE) +
  geom_richtext(data=tibble(x=5.5, y=22), fill = NA, label.color = NA, label="r = 0.65",aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=5.5, y=21.5), fill = NA, label.color = NA, label="*p* = 0.00",aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  labs(x = "Log(CFU/g lung)", y = bquote("Log(CCL2/gram lung)")) +
  scale_shape_manual(breaks = c("c57", "het", "dbdb"),
                     values = c(22,21,21),
                     labels = c("+/+", "db/+", "db/db")) +
  scale_fill_manual(breaks = c("c57", "het", "dbdb"),
                    values = c("black","#5F4B8B","#AD5E99"),
                    labels = c("+/+", "db/+", "db/db")) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10))
ggsave("graphs/lung_CCL2_corr.pdf", height = 6, width = 6, unit = "cm")
#TNF
rcorr(TNF_exp_stat$log_CFU_per_g_total, TNF_exp_stat$log_TNF_pg_per_g)
ggplot(TNF_exp_stat, aes(x = log_TNF_pg_per_g, y = log_CFU_per_g_total, fill = Genotype, shape = Genotype)) +
  geom_point(size = 2.5,show.legend = FALSE) +
  geom_richtext(data=tibble(x=5, y=22), fill = NA, label.color = NA, label="r = 0.4",aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=5, y=21.5), fill = NA, label.color = NA, label="*p* = 0.01",aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  labs(x = "Log(CFU/g lung)", y = bquote("Log("*TNF*alpha*"/gram lung)")) +
  scale_shape_manual(breaks = c("c57", "het", "dbdb"),
                     values = c(22,21,21),
                     labels = c("+/+", "db/+", "db/db")) +
  scale_fill_manual(breaks = c("c57", "het", "dbdb"),
                    values = c("black","#5F4B8B","#AD5E99"),
                    labels = c("+/+", "db/+", "db/db")) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10))
ggsave("graphs/lung_TNF_corr.pdf", height = 6, width = 6, unit = "cm")
#VEGF
rcorr(VEGF_exp_stat$log_CFU_per_g_total, VEGF_exp_stat$log_VEGF_pg_per_g)
ggplot(VEGF_exp_stat, aes(x = log_VEGF_pg_per_g, y = log_CFU_per_g_total, fill = Genotype, shape = Genotype)) +
  geom_point(size = 2.5,show.legend = FALSE) +
  geom_richtext(data=tibble(x=4.25, y=22), fill = NA, label.color = NA, label="r = 0.28",aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  geom_richtext(data=tibble(x=4.25, y=21.5), fill = NA, label.color = NA, label="*p* = 0.09",aes(x=x, y=y), inherit.aes=FALSE, size=2.5) +
  labs(x = "Log(CFU/g lung)", y = bquote("Log(VEGF-A/gram lung)")) +
  scale_shape_manual(breaks = c("c57", "het", "dbdb"),
                     values = c(22,21,21),
                     labels = c("+/+", "db/+", "db/db")) +
  scale_fill_manual(breaks = c("c57", "het", "dbdb"),
                    values = c("black","#5F4B8B","#AD5E99"),
                    labels = c("+/+", "db/+", "db/db")) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_markdown(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "lightgrey"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "lightgrey"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10))
ggsave("graphs/lung_VEGF_corr.pdf", height = 6, width = 6, unit = "cm")

