# preservation  module comparisons 

library(tidyverse)

datExpr <- readRDS("manuscript/brain/manuscript70/results/WGCNA/datExpr_MEA_CDOM.RDS") 
net <- readRDS("manuscript/brain/manuscript70/results/WGCNA/net_MEA_CDOM_Power9.RDS")
MEs <- readRDS("manuscript/brain/manuscript70/results/wgcna/MEA_CDOM_MEs.RDS")
coldata <- read_csv("manuscript/brain/sample70min_table.csv")
head(coldata)
str(coldata)

#group thing 
coldata$groupEX <- coldata$group
# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))

coldata <-coldata %>% mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))
coldata <- coldata %>%  dplyr::select(-group, -period)

#getting condition1
table(coldata$condition)
# CDOM, RDOM to Descenders (DOM to SUB)(4->1)
# CSUB, SUB to Ascenders (Sub to DOM)  (1->4)

coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "DOM", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "SUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "DES", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "ASC", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 4, "CSUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 1, "CDOM", coldata$condition1)

coldata <- coldata %>% 
  filter(Postrank != 3) %>% 
  filter(condition1 != 'ascenders') %>% 
  filter(region == "AMY")



ME_dom <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleName") %>%
  pivot_longer(cols = 2:23, names_to = "Module") %>% 
  full_join(coldata)

head(ME_dom)

ME_dom$condition1

ME_dom$Module <- gsub("ME", "", ME_dom$Module)


ME_dom  <- ME_dom %>%
  mutate(status = condition1) %>% 
  filter(condition1 != "ASC") %>%
  filter(condition1 != "SUB") %>%
  filter(condition1 != "CSUB")

ME_dom$status <- factor(ME_dom$status, levels = c("CDOM", "DOM", "DES"))

#SUBs 
datExpr <- readRDS("manuscript/brain/manuscript70/results/WGCNA/datExpr_MEA_CSUB.RDS") 
net <- readRDS("manuscript/brain/manuscript70/results/WGCNA/net_MEA_CSUB_Power10.RDS")
MEs <- readRDS("manuscript/brain/manuscript70/results/wgcna/MEA_CSUB_MEs.RDS")

ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleName") %>%
  pivot_longer(cols = 2:24, names_to = "Module") %>% 
  full_join(coldata)

head(ME_df)

ME_df$condition1
ME_df <- ME_df %>%
  mutate(status = condition1) %>% 
  filter(condition1 != "DES") %>%
  filter(condition1 != "DOM") %>%
  filter(condition1 != "CDOM")  %>% 
  filter(condition1 != "ascenders") %>%
  filter(condition1 != "descenders") %>%
  filter(condition1 != "control") %>%
  filter(condition1 != "same") 
ME_df$Module <- gsub("ME", "", ME_df$Module)

ME_df$status <- factor(ME_df$status, levels = c("CSUB", "SUB", "ASC"))

### comparison boxplots for mPFC 

source("functions/geom_boxjitter.R")

# red <- ME_dom %>% filter(Module == "red")
brown <-ME_dom %>% filter(Module == "brown") 
# yellow <-ME_dom %>% filter(Module == "yellow") 
pink <-ME_dom %>% filter(Module == "pink") 
green<-ME_dom %>% filter(Module == "green")
blue <- ME_dom %>% filter(Module == "blue")
purple <- ME_dom %>% filter(Module == "purple")



turq <- ME_df %>% filter(Module == "turquoise")
red <- ME_df %>% filter(Module == "red")
subg <- ME_df %>% filter(Module == "green")
suby <- ME_df %>% filter(Module == "yellow")
sal <- ME_df %>% filter(Module == "salmon")
red %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Red Module")+ 
  theme_classic()+
  theme(legend.position = "none",  text=element_text(size=15))  -> p_red


subg %>% 
ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Green Module")+ 
  theme_classic()+
  theme(legend.position = "none",  text=element_text(size=15))  -> p_subg

brown %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Brown Module")+ 
  theme_classic()+
  theme(legend.position = "none", text=element_text(size=15))  -> p_brown

sal %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Salmon Module")+ 
  theme_classic()+
  theme(legend.position = "none", text=element_text(size=15))  -> p_sal

yellow %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Yellow Module")+ 
  theme_classic()+
  theme(legend.position = "none",  text=element_text(size=15))  -> p_yellow

green %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Green Module")+ 
  theme_classic()+
  theme(legend.position = "none", text=element_text(size=15))  -> p_green

blue %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Blue Module")+ 
  theme_classic()+
  theme(legend.position = "none",  text=element_text(size=15))  -> p_blue




purple %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Purple Module")+ 
  theme_classic()+
  theme(legend.position = "none",  text=element_text(size=15))  -> p_purple




turq %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Turquoise Module")+ 
  theme_classic() + theme(legend.position = "none",  text=element_text(size=15))  -> p_turq



suby %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Yellow Module")+ 
  theme_classic()+
  theme(legend.position = "none", text=element_text(size=15))  -> p_suby

lg %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "LightGreen Module")+ 
  theme_classic()+
  theme(legend.position = "none", text=element_text(size=15))  -> p_lg

pink %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Pink Module")+ 
  theme_classic()+
  theme(legend.position = "none", text=element_text(size=15))  -> p_pink


cp <- gridExtra::grid.arrange(p_pink, p_sal, p_purple, p_turq, p_brown, p_suby,p_green, p_subg, nrow =4) 

ggsave("manuscript/brain/manuscri|pt70/results/results_figures/mPFC_interesting_preservations.png", cp, 
       width = 8, height = 8,dpi = 150)

# others
noncp <- gridExtra::grid.arrange(p_brown, p_turq,p_pink, p_green,p_lg, p_suby, nrow =3)


ggsave("manuscript/brain/manuscript70/results/results_figures/mPFC_preservations_NOTusing.png", noncp, 
       width = 12, height = 12,dpi = 150)

### checking significants 
lm_dom <- readRDS("manuscript/brain/manuscript70/results/RDS/mPFC_CDOM_lm_result_all.RDS")
lm_sub <- readRDS("manuscript/brain/manuscript70/results/RDS/mPFC_CSUB_lm_result_all.RDS")


#########MEA
library(tidyverse)
library(WGCNA)
datExpr <- readRDS("manuscript/brain/manuscript70/results/WGCNA/datExpr_MEA_CDOM.RDS") 
net <- readRDS("manuscript/brain/manuscript70/results/WGCNA/net_MEA_CDOM_Power9.RDS")
MEs <- readRDS("manuscript/brain/manuscript70/results/wgcna/MEA_CDOM_MEs.RDS")


coldata <- read_csv("manuscript/brain/sample70min_table.csv")
head(coldata)
str(coldata)


coldata$groupEX <- coldata$group
# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))

coldata <-coldata %>% mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))
coldata <- coldata %>%  dplyr::select(-group, -period)

#getting condition1
table(coldata$condition)
# CDOM, RDOM to Descenders (DOM to SUB)(4->1)
# CSUB, SUB to Ascenders (Sub to DOM)  (1->4)

coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "DOM", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "SUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "DES", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "ASC", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 4, "CSUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 1, "CDOM", coldata$condition1)

coldata <- coldata %>% 
  filter(Postrank != 3) %>% 
  filter(condition1 != 'ascenders') %>% 
  filter(region == "AMY")



ME_dom <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleName") %>%
  pivot_longer(cols = 2:23, names_to = "Module") %>% 
  full_join(coldata)

head(ME_dom)

ME_dom$condition1

ME_dom$Module <- gsub("ME", "", ME_dom$Module)


ME_dom  <- ME_dom %>%
  mutate(status = condition1) %>% 
  filter(condition1 != "ASC") %>%
  filter(condition1 != "SUB") %>%
  filter(condition1 != "CSUB")

ME_dom$status <- factor(ME_dom$status, levels = c("CDOM", "DOM", "DES"))


#SUBs 
datExpr <- readRDS("manuscript/brain/manuscript70/results/WGCNA/datExpr_MEA_CSUB.RDS") 
net <- readRDS("manuscript/brain/manuscript70/results/WGCNA/net_MEA_CSUB_Power10.RDS")
MEs <- readRDS("manuscript/brain/manuscript70/results/wgcna/MEA_CSUB_MEs.RDS")


ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleName") %>%
  pivot_longer(cols = 2:25, names_to = "Module") %>% 
  full_join(coldata)

head(ME_df)

ME_df$condition1
ME_df <- ME_df %>%
  mutate(status = condition1) %>% 
  filter(condition1 != "DES") %>%
  filter(condition1 != "DOM") %>%
  filter(condition1 != "CDOM")  %>% 
  filter(condition1 != "ascenders") %>%
  filter(condition1 != "descenders") %>%
  filter(condition1 != "control") %>%
  filter(condition1 != "same") 
ME_df$Module <- gsub("ME", "", ME_df$Module)

ME_df$status <- factor(ME_df$status, levels = c("CSUB", "SUB", "ASC"))



brown <-ME_dom %>% filter(Module == "brown") 
purple <-ME_dom %>% filter(Module == "purple") 
pink <-ME_dom %>% filter(Module == "pink") 
blue<-ME_dom %>% filter(Module == "blue")
lg <-ME_dom %>% filter(Module == "greenyellow") 
green<-ME_dom %>% filter(Module == "green")

turq <- ME_df %>% filter(Module == "turquoise")
red <- ME_df %>% filter(Module == "red")
salmon <- ME_df %>% filter(Module == "salmon")
suby <- ME_df %>% filter(Module == "yellow")
subg <- ME_df %>% filter(Module == "green")


turq %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Turquoise Module")+ 
  theme_classic() + theme(legend.position = "none")  -> p_turq



red %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Red Module")+ 
  theme_classic()+
  theme(legend.position = "none")  -> p_red


salmon %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Salmon Module")+ 
  theme_classic() + theme(legend.position = "none")  -> p_salmon



suby %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Yellow Module")+ 
  theme_classic()+
  theme(legend.position = "none")  -> p_suby

subg %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Green Module")+ 
  theme_classic()+
  theme(legend.position = "none")  -> p_subg



purple %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Purple Module")+ 
  theme_classic()+
  theme(legend.position = "none")  -> p_purple

pink %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Pink Module")+ 
  theme_classic()+
  theme(legend.position = "none")  -> p_pink


brown %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Brown Module")+ 
  theme_classic()+
  theme(legend.position = "none")  -> p_brown

blue %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Blue Module")+ 
  theme_classic()+
  theme(legend.position = "none")  -> p_blue

green %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Green Module")+ 
  theme_classic()+
  theme(legend.position = "none")  -> p_green
lg %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "GreenYellow Module")+ 
  theme_classic()+
  theme(legend.position = "none")  -> p_gl


gridExtra::grid.arrange(p_brown, p_suby, nrow =1)# interesting
gridExtra::grid.arrange(p_pink, p_red, nrow =1) #?
gridExtra::grid.arrange(p_pink, p_salmon, nrow =1)# interesting
gridExtra::grid.arrange(p_blue, p_turq, nrow =1)# ?
gridExtra::grid.arrange(p_purple, p_turq, nrow =1)# interesting
gridExtra::grid.arrange(p_green, p_subg, nrow =1) # interesting
gridExtra::grid.arrange(p_gl, p_subg, nrow =1)

cp <- gridExtra::grid.arrange(p_pink, p_salmon,p_purple, p_turq, p_brown, p_suby,p_green, p_subg, nrow =4) 

ggsave("manuscript/brain/manuscript70/results/results_figures/MEA_interesting_preservations.png", cp, 
       width = 10, height = 14,dpi = 150)

# others
noncp <- gridExtra::grid.arrange(p_pink, p_red,p_blue, p_turq,p_gl, p_subg, nrow =3)


ggsave("manuscript/brain/manuscript70/results/results_figures/MEA_preservations_NOTusing.png", noncp, 
       width = 12, height = 12,dpi = 150)

### checking significants 
lm_dom <- readRDS("manuscript/brain/manuscript70/results/RDS/MEA_CDOM_lm_result_all.RDS")
lm_sub <- readRDS("manuscript/brain/manuscript70/results/RDS/MEA_CSUB_lm_result_all.RDS")
