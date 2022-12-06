library(limma)
library(Mus.musculus)
library(DESeq2)
library(tidyverse)
global_size = 10


#Getting metadata ready 
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

amy_data <- coldata %>% 
  filter(Postrank != 3) %>% 
  filter(condition1 != 'ascenders') %>% 
  filter(condition1 != "ASC") %>% 
  filter(condition1 != "CSUB")%>% 
  filter(condition1 != "SUB")%>% 
  filter(region == "AMY")

#just get samples I want
row.names <- amy_data$SampleName

amy_data1 <- amy_data%>% dplyr::select(-SampleName)

row.names(amy_data1) <- row.names #Assigning row names from as sample names  
head(amy_data1)

# Bring in count data for mAMY
a_countdata <- read_csv("manuscript/brain/AMY_counts.csv")
a_countdata$X -> nrows
a_count <- a_countdata[,-1]
row.names(a_count) <- nrows
a_count <-as.data.frame(a_count)

#filter
a_count <- a_count[rowSums(a_count >= 25), ]

#checks
all (row.names(amy_data1) %in% colnames(a_count)) #check 
a_count <- a_count[,rownames(amy_data1)]
all(rownames(amy_data1) == colnames(a_count)) #check



# Second mpfc data data 
pfc_data<- coldata %>% filter(region != "AMY")%>% 
  filter(condition1 != "ASC") %>% 
  filter(condition1 != "CSUB")%>% 
  filter(condition1 != "SUB")%>% 
  filter(condition1 != "same")%>% 
  filter(condition1 != "ascenders") %>% 
  filter(condition1 != "descenders")%>%
  filter(condition1 != "control")%>%
  arrange(condition1) 

row.names <- pfc_data$SampleName

pfc_data1 <- pfc_data%>% dplyr::select(-SampleName)

row.names(pfc_data1) <- row.names #Assigning row names from as sample names  
head(pfc_data1)

# Bring in count data for mpfc
p_countdata <- read_csv("manuscript/brain/PFC_counts.csv")
p_countdata$X -> nrows
p_count <- p_countdata[,-1]
row.names(p_count) <- nrows
p_count <-as.data.frame(p_count)

#filter
p_count <- p_count[rowSums(p_count >= 25), ]
# pkeep <- rowSums(p_count >= 10)
# p_count <- p_count[pkeep,]


p_count <- p_count[, rownames(pfc_data1)]
all(rownames(pfc_data1) == colnames(p_count)) #check


#MeA
dag <- DESeqDataSetFromMatrix(countData = a_count,
                              colData = amy_data1,
                              design = ~condition1)



nrow(dag) # 11951 lost about half

dag_sub <- dag[, dag$condition1 %in% c("CDOM", "DOM", "DES")]
dag_sub$condition1 <- droplevels(dag_sub$condition1)
dag_sub$condition1
dag<- dag_sub
dag$condition1



#mPFC
dpg <- DESeqDataSetFromMatrix(countData = p_count,
                              colData = pfc_data1,
                              design = ~condition1)

nrow(dpg) # 12051
dpg_sub <- dpg[, dpg$condition1 %in% c("CDOM", "DOM", "DES")]
dpg_sub$condition1 <- droplevels(dpg_sub$condition1)
dpg_sub$condition1
dpg<- dpg_sub
dpg$condition1


#dds
da_dom <- DESeq(dag)
dp_dom <- DESeq(dpg)



vsd <- vst(da_dom, blind=FALSE) 
plotPCA(vsd, intgroup = c("condition1"), returnData = T, ntop=200) %>%
  left_join(amy_data %>% 
              rename(name = SampleName))-> d

d$condition1 <- factor(d$condition1, levels = c("CDOM", "DOM", "DES"))
d$group <- factor(d$group, levels = c("CDOM", "DOM", "DES"))
ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
  geom_point(size = 3, alpha = 0.5) + 
  # geom_text(aes(label = subjectID), vjust = -0.5)+
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  labs(x = paste0("PC1: ", 
                  round(attributes(d)$percentVar[1] * 100), 
                  "% variance"),
       y = paste0("PC2: ", 
                  round(attributes(d)$percentVar[2] * 100), 
                  "% variance"),
       color = "Status",
       fill = "Status") +
  theme_classic()+
  theme(legend.position = c(0.15,0.20),
  legend.key.size = unit(0.3, 'cm'),
  legend.text = element_text(size=12),
  legend.title = element_text(size=12),
  text=element_text(size=15))-> p

p

dev.off()
ggsave("manuscript/brain/manuscript70/results/results_figures/CDOM_MEA70_PCA_man.png",p,width = 4.75, height = 4.25, dpi=300)




library(viridis)
library(scales)

nCol <- 50
myCol <- viridis(n = nCol)
myCol


vsd <- vst(dp_dom, blind=FALSE)
plotPCA(vsd, intgroup = c("condition1"), returnData = T,ntop=200) %>%
  left_join(pfc_data %>% 
              rename(name = SampleName))-> d
d$group <- factor(d$group, levels =c("CDOM", "DOM", "DES"))  

d2 <- d %>% filter(post_idbatch != "3-4Batch12")

ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
  geom_point(size = 3, alpha = 0.5) + 
  # geom_text(aes(label =), vjust = -0.5)+
  scale_color_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  labs(x = paste0("PC1: ", 
                  round(attributes(d)$percentVar[1] * 100), 
                  "% variance"),
       y = paste0("PC2: ", 
                  round(attributes(d)$percentVar[2] * 100), 
                  "% variance"),
       color = "Status",
       fill = "Status") +
  theme_classic()+
  theme(legend.position = c(0.15,0.85),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        text=element_text(size=15)) -> p1

p
p1
dev.off()

ggsave("manuscript/brain/manuscript70/results/results_figures/CDOM_mPFC70_PCA_man.png",p1,width = 4.75, height = 4.25, dpi=300)



#Getting metadata ready 
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

amy_data <- coldata %>% 
  filter(Postrank != 3) %>% 
  filter(condition1 != 'ascenders') %>% 
  filter(condition1 != "DOM") %>% 
  filter(condition1 != "CDOM")%>% 
  filter(condition1 != "DES")%>% 
  filter(region == "AMY")

#just get samples I want
row.names <- amy_data$SampleName

amy_data1 <- amy_data%>% dplyr::select(-SampleName)

row.names(amy_data1) <- row.names #Assigning row names from as sample names  
head(amy_data1)

# Bring in count data for mAMY
a_countdata <- read_csv("manuscript/brain/AMY_counts.csv")
a_countdata$X -> nrows
a_count <- a_countdata[,-1]
row.names(a_count) <- nrows
a_count <-as.data.frame(a_count)

#filter
a_count <- a_count[rowSums(a_count >= 25), ]

#checks
all (row.names(amy_data1) %in% colnames(a_count)) #check 
a_count <- a_count[,rownames(amy_data1)]
all(rownames(amy_data1) == colnames(a_count)) #check



# Second mpfc data data 
pfc_data<- coldata %>% filter(region != "AMY")%>% 
  filter(condition1 != "DES") %>% 
  filter(condition1 != "CDOM")%>% 
  filter(condition1 != "DOM")%>% 
  filter(condition1 != "same")%>% 
  filter(condition1 != "ascenders") %>% 
  filter(condition1 != "descenders")%>%
  filter(condition1 != "control")%>%
  arrange(condition1) 

row.names <- pfc_data$SampleName

pfc_data1 <- pfc_data%>% dplyr::select(-SampleName)

row.names(pfc_data1) <- row.names #Assigning row names from as sample names  
head(pfc_data1)

# Bring in count data for mpfc
p_countdata <- read_csv("manuscript/brain/PFC_counts.csv")
p_countdata$X -> nrows
p_count <- p_countdata[,-1]
row.names(p_count) <- nrows
p_count <-as.data.frame(p_count)

#filter
p_count <- p_count[rowSums(p_count >= 25), ]
# pkeep <- rowSums(p_count >= 10)
# p_count <- p_count[pkeep,]


p_count <- p_count[, rownames(pfc_data1)]
all(rownames(pfc_data1) == colnames(p_count)) #check


#MeA
dag <- DESeqDataSetFromMatrix(countData = a_count,
                              colData = amy_data1,
                              design = ~condition1)



nrow(dag) # 11951 lost about half

dag_sub <- dag[, dag$condition1 %in% c("CSUB", "SUB", "ASC")]
dag_sub$condition1 <- droplevels(dag_sub$condition1)
dag_sub$condition1
dag<- dag_sub
dag$condition1



#mPFC
dpg <- DESeqDataSetFromMatrix(countData = p_count,
                              colData = pfc_data1,
                              design = ~condition1)

nrow(dpg) # 12051
dpg_sub <- dpg[, dpg$condition1 %in% c("CSUB", "SUB", "ASC")]
dpg_sub$condition1 <- droplevels(dpg_sub$condition1)
dpg_sub$condition1
dpg<- dpg_sub
dpg$condition1


#dds
da_dom <- DESeq(dag)


dp_dom <- DESeq(dpg)


vsd <- vst(da_dom, blind=FALSE) 
plotPCA(vsd, intgroup = c("condition1"), returnData = T, ntop=200) %>%
  left_join(amy_data %>% 
              rename(name = SampleName))-> d

d$group <- factor(d$group, levels = c("CSUB", "SUB", "ASC"))
ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
  geom_point(size = 3, alpha = 0.5) + 
  # geom_text(aes(label = subjectID), vjust = -0.5)+
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  labs(x = paste0("PC1: ", 
                  round(attributes(d)$percentVar[1] * 100), 
                  "% variance"),
       y = paste0("PC2: ", 
                  round(attributes(d)$percentVar[2] * 100), 
                  "% variance"),
       color = "Status",
       fill = "Status") +
  theme_classic()+
  theme(legend.position = c(0.15,0.20),
  legend.key.size = unit(0.4, 'cm'),
  legend.text = element_text(size=12),
  legend.title = element_text(size=12),
  text=element_text(size=15)) -> p

p
dev.off()
ggsave("manuscript/brain/manuscript70/results/results_figures/CSUB_MEA70_PCA_man.png",p,width = 4.75, height = 4.25, dpi=300)

vsd <- vst(dp_dom, blind=FALSE)
plotPCA(vsd, intgroup = c("condition1"), returnData = T,ntop=200) %>%
  left_join(pfc_data %>% 
              rename(name = SampleName))-> d
d$group <- factor(d$group, levels =c("CSUB", "SUB", "ASC"))  


ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
  geom_point(size = 3, alpha = 0.5) + 
  # geom_text(aes(label =), vjust = -0.5)+
  scale_color_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  labs(x = paste0("PC1: ", 
                  round(attributes(d)$percentVar[1] * 100), 
                  "% variance"),
       y = paste0("PC2: ", 
                  round(attributes(d)$percentVar[2] * 100), 
                  "% variance"),
       color = "Status",
       fill = "Status") +
  theme_classic()+
  theme(legend.position = c(0.15,0.85),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        text=element_text(size=15)) -> p1


p
p1

dev.off()
ggsave("manuscript/brain/manuscript70/results/results_figures/CSUB_mPFC70_PCA_man.png",p1,width = 4.75, height = 4.25, dpi=300)


# venn Diagram 
library(vennplot)
library(systemPipeR)

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_mPFC_ControlDD.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= .2)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 




deg.up <- limma_list %>% map(~filter(.,logFC >= 0.2)) %>% map(~dplyr::select(.,symbol))
# up.cort <- mPFC_cort %>% filter(.,logFC >= 0.2) %>% dplyr::select(.,symbol)
l.up <- list("CON-DOM" = c(t(deg.up$controldom)),"CON-DES" = c(t(deg.up$controldes)),"DOM-DES"= c(t(deg.up$domdes)))
vennset.up <- overLapper(l.up[1:3], type="vennsets")
vennPlot(vennset.up)

deg.down <- limma_list %>% map(~filter(.,logFC <= 0.2)) %>% map(~dplyr::select(.,symbol))
# down.cort <- mPFC_cort %>% filter(.,logFC <= 0.2) %>% dplyr::select(.,symbol)
l.down <- list("CON-DOM" = c(t(deg.down$controldom)),"CON-DES" = c(t(deg.down$controldes)),"DOM-DES"= c(t(deg.down$domdes)))
vennset.down <- overLapper(l.down[1:3], type="vennsets")
vennPlot(vennset.down)

vennPlot(list(vennset.down,vennset.up),mymain="", mysub="", colmode=2, ccol=c( "red", "blue"))



##subs

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_MeA_ControlSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= .2)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 


# venn Diagram 
library(vennplot)
library(systemPipeR)
deg.up <- limma_list %>% map(~filter(.,logFC >= 0.2)) %>% map(~dplyr::select(.,symbol))
# up.cort <- mPFC_cort %>% filter(.,logFC >= 0.2) %>% dplyr::select(.,symbol)
l.up <- list("CON-SUB" = c(t(deg.up$controlsub)),"CON-ASC" = c(t(deg.up$controlasc)),"SUB-ASC"= c(t(deg.up$subasc)))
vennset.up <- overLapper(l.up[1:3], type="vennsets")
vennPlot(vennset.up)

deg.down <- limma_list %>% map(~filter(.,logFC <= 0.2)) %>% map(~dplyr::select(.,symbol))
# down.cort <- mPFC_cort %>% filter(.,logFC <= 0.2) %>% dplyr::select(.,symbol)
l.down <- list("CON-SUB" = c(t(deg.down$controlsub)),"CON-ASC" = c(t(deg.down$controlasc)),"SUB-ASC"= c(t(deg.down$subasc)))
vennset.down <- overLapper(l.down[1:3], type="vennsets")
vennPlot(vennset.down)

vennPlot(list(vennset.down,vennset.up),mymain="", mysub="", colmode=2, ccol=c( "red", "blue"))





