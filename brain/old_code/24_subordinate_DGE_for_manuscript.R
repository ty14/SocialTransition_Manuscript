library(Mus.musculus)
organism = 'org.Mm.eg.db'
library(organism, character.only = TRUE)
library(biomaRt)
library(AnnotationDbi)
library(pheatmap)
library(annotables)
library(clusterProfiler)
library(enrichplot)
library(organism, character.only = TRUE)
library(DOSE)
library(VennDiagram)
library(limma)
library(tidyverse)
grcm38 # mouse genes
library("ggsignif")  


my_logFC_threshold = 0.2

#MEA data

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_MEA_ControlSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

cs <- limma_list$controlsub
cs$reg <- ifelse(cs$logFC >= 0.2, "up", "down")

ca <- limma_list$controlasc 
ca$reg <- ifelse(ca$logFC >= 0.2, "up", "down")

sa <- limma_list$subasc 
sa$reg <- ifelse(sa$logFC >= 0.2, "up", "down")

#looking for overlap: 
cdd <- cs %>% filter(reg == "down")

ddd <- sa %>% filter(reg == "up")

up_sub <- cdd$symbol[(cdd$symbol %in% ddd$symbol)] %>% as.data.frame(.)
colnames(up_sub)<-"symbol"

up_sub<- up_sub %>% mutate(reg = "Up")

cddd <- cs %>% filter(reg == "up")

dddd <- sa %>% filter(reg == "down")

down_sub <- cddd$symbol[(cddd$symbol %in% dddd$symbol)] %>% as.data.frame(.)
colnames(down_sub)<-"symbol"

down_sub <- down_sub %>% mutate(reg = "Down")

#up dom = 151, down dom = 124
sub <- up_sub%>% full_join(down_sub)



write.csv(sub, 'manuscript/brain/manuscript70/results/tables/subordinate_distruption_MEA_genes.csv' ,row.names = F)

#looking further into staying subordinate  genes:
sub<- read_csv('manuscript/brain/manuscript70/results/tables/subordinate_distruption_MEA_genes.csv')
head(sub)


casc_top <- cs[cs$symbol %in% sub$symbol,]

cd_up <- casc_top %>% arrange(logFC)
x <- head(cd_up, 15)

sa_top <- sa[sa$symbol %in% sub$symbol,]
sa_up <- sa_top %>% arrange(-logFC)
y <- head(sa_up, 15)


x[x$symbol %in% y$symbol,]

#UP:Sim1, Neurod6, Cbln2, Rgs13,  Dach2, Tmem40

#DOWN:  Osr1, Pifo, Gfap, Alox12b,  Tap1, Htr4 


### Go terms 
my_ont = "BP"
my_showCategory = 100

gettop10GO <- function(limma_df,my_showCategory = 10){
  go_df_up <- limma_df %>% filter(reg == "Up") %>% 
    left_join(grcm38 %>% dplyr::select(symbol, entrez)) %>% 
    mutate(entrez = as.character(entrez)) %>% 
    dplyr::select(entrez) %>% 
    filter(!is.na(entrez))
  
  
  go_df_down <- limma_df %>% filter(reg == "Down") %>% 
    left_join(grcm38 %>% dplyr::select(symbol, entrez)) %>% 
    mutate(entrez = as.character(entrez)) %>% 
    dplyr::select(entrez) %>% 
    filter(!is.na(entrez)) 
  
  ggo_up <- enrichGO(gene = go_df_up$entrez %>% unique(),
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                     keyType = "ENTREZID",
                     ont = my_ont,
                     readable = T,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.2,
                     qvalueCutoff  = 0.20)
  
  ggo_down <- enrichGO(gene = go_df_down$entrez %>% unique() ,
                       OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = my_ont,
                       readable = T,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.2, 
                       qvalueCutoff  = 0.20) 
  
  fortify(
    ggo_up,
    showCategory = my_showCategory,
    by = "GeneRatio",
    split = NULL,
    includeAll = TRUE) %>% 
    dplyr::arrange(desc(GeneRatio)) %>% 
    dplyr::mutate(direction = "Up") -> temp1
  
  
  fortify(
    ggo_down,
    showCategory = my_showCategory,
    by = "GeneRatio",
    split = NULL,
    includeAll = TRUE) %>% 
    dplyr::arrange(desc(GeneRatio)) %>% 
    dplyr::mutate(direction = "Down") -> temp2
  return(rbind(temp1,temp2))
  
}


gettop10GO(sub, my_showCategory) -> top10go1

write.csv(top10go1, 'manuscript/brain/manuscript70/results/tables/SUB_GOterms_BP.csv' ,row.names = F)


top10go1$Description[1:5]
top10go1$Description[6:10]
top10go1$Description[11:15]
top10go1$Description[16:20]
top10go1$Description[21:25]

up <- top10go1 %>% arrange(qvalue) %>% filter(direction == "Up")
head(up$Description)


down <- top10go1 %>% arrange(qvalue) %>% filter(direction == "Down")
head(down$Description)



### Make DGE boxplots:

ex <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_MeA_ControlSUB.RDS")

sm <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

sub_m_id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

# Sim1, Neuro6, Cbln2, Rgs13,  Dach2, Tmem40

a1<- sm %>% filter(symbol == "Sim1")
a2 <- sm %>% filter(symbol == "Neurod6")
a3 <- sm %>% filter(symbol == "Cbln2")
a4 <- sm %>% filter(symbol == "Rgs13")
a5<- sm %>% filter(symbol == "Dach2")
a6 <- sm %>% filter(symbol == "Tmem40")


dp <- a1 %>% rbind(a2,a3,a4,a5,a6)


colnames(dp)

dp2 <- dp  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(sub_m_id)

dp2$group <- factor(dp2$group, levels = c("CSUB", "SUB", "ASC")) 


source('functions/geom_boxjitter.R')


dp.l <- split(dp2 , dp2$symbol)
#Sim1 c *, d **
# Neurod6,c*, d **
#Cbln2, c ***, d ***
#Rgs13,c **, d **
#Dach2 c **, **
#Tmem40 c **, d **
# MEA Normalized Expression

c1 <- dp.l[[5]]
c1p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("MEA Normalized Expression") +
  xlab("")+
  ggtitle(substitute(paste(italic("Sim1"))))+
  ylim (-3,9.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(8.5, 7.5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "**"))

c1p


c1 <- dp.l[[3]]
c2p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Neurod6"))))+
  ylim ( 2,8)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(7.5, 6.6),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "**"))

c2p





c1 <- dp.l[[1]]
c3p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Cbln2"))))+
  ylim ( 2,10.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(9.5,8.6),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("***", "***"))

c3p


c1 <- dp.l[[4]]
c4p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Rgs13"))))+
  ylim ( -4,6)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(5, 4),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "**"))

c4p


c1 <- dp.l[[2]]
c5p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Dach2"))))+
  ylim ( 1,6.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(6, 5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "**"))

c5p


c1 <- dp.l[[6]]
c6p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Tmem40"))))+
  ylim ( -3, 5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(4, 3),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "**"))

c6p


subup <- gridExtra::grid.arrange(c1p, c2p, c5p, c4p, c3p, c6p, nrow = 1)

ggsave("manuscript/brain/manuscript70/results/results_figures/SUB_TOPRegulated_genes_manuscriptSignifcants.png",domup, height = 3, width= 18,  dpi=300)

# Osr1, Pifo, Gfap, Alox12b,  Tap1, Htr4 

a1<- sm %>% filter(symbol == "Osr1")
a2 <- sm %>% filter(symbol == "Pifo")
a3 <- sm %>% filter(symbol == "Gfap")
a4 <- sm %>% filter(symbol == "Alox12b")
a5<- sm %>% filter(symbol == "Tap1")
a6 <- sm %>% filter(symbol == "Htr4")


dp <- a1 %>% rbind(a2,a3,a4,a5,a6)


colnames(dp)

dp2 <- dp  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(sub_m_id)

dp2$group <- factor(dp2$group, levels = c("CSUB", "SUB", "ASC")) 


source('functions/geom_boxjitter.R')


dp.l <- split(dp2 , dp2$symbol)
#Osr1, *, *
# Pifo, *, * 
# Gfap, *. *
# Alox12b, *, * 
# Tap1, *, *
# Htr4 **, ***
# MEA Normalized Expression

c1 <- dp.l[[5]]
c1p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("MEA Normalized Expression") +
  xlab("")+
  ggtitle(substitute(paste(italic("Pifo"))))+
  ylim (-3,5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(4, 3),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))

c1p


c1 <- dp.l[[3]]
c2p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Htr4"))))+
  ylim ( -1,6)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(5, 4),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "***"))

c2p





c1 <- dp.l[[1]]
c3p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Alox12b"))))+
  ylim ( -3,5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(4.25, 3.25),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))

c3p


c1 <- dp.l[[4]]
c4p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Osr1"))))+
  ylim ( -3,6)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(5.25, 4.25),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))

c4p


c1 <- dp.l[[2]]
c5p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Gfap"))))+
  ylim ( 3.5,10)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(9.5, 8.5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))

c5p


c1 <- dp.l[[6]]
c6p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Tap1"))))+
  ylim ( -1.5, 5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(4.25, 3.25),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))

c6p


subdown <- gridExtra::grid.arrange(c1p, c2p, c5p, c4p, c3p, c6p, nrow = 1)

ggsave("manuscript/brain/manuscript70/results/results_figures/SUB_DOWN_TOPRegulated_genes_manuscriptSignifcants.png",subdown, height = 3, width= 18,  dpi=300)



########## cluster go terms stuff 

df <- read_csv('manuscript/brain/manuscript70/results/tables/SUB_GOterms_BP.csv')

head(df) 

#up data 
df_up <- df %>%  filter(direction == "Up") 
df_upx <- df_up[1:30,]

# down data 
df_down <- df %>%  filter(direction == "Down") 
df_downx <- df_down[1:30,]

#get columns we want - check count to make sure it is greater than 5. 
df_upx <-  df_upx[c("ID", "Description", "pvalue", "qvalue", "geneID")]

df_downx <-   df_downx[c("ID", "Description", "pvalue", "qvalue", "geneID")]

# format gene list column
df_upx$geneID <- gsub("/", ",",  df_upx$geneID)

df_downx$geneID <- gsub("/", ",",  df_downx$geneID)

# add column for phenotype
df_upx <- cbind( df_upx, phenotype= 1)
df_upx <-  df_upx[, c(1, 2, 3, 4, 6, 5)]

df_downx <- cbind( df_downx, phenotype= -1)
df_downx <-  df_downx[, c(1, 2, 3, 4, 6, 5)]
# change column headers
colnames( df_upx) <- c("Name","Description", "pvalue","qvalue", "phenotype", "genes")
head(df_upx)
# write file for cytoscape
df_upx.filename <-file.path("manuscript/brain/manuscript70/results/tables/MEA_SUB_up30_clusterprofiler_cluster_enr_results.txt")
write.table(df_upx,df_upx.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)



colnames(df_downx) <- c("Name","Description", "pvalue","qvalue", "phenotype", "genes")
head(df_downx)
# write file for cytoscape
df_downx.filename <-file.path("manuscript/brain/manuscript70/results/tables/MEA_SUB_down30_clusterprofiler_cluster_enr_results.txt")
write.table(df_downx,df_downx.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)



########### heat plot

sub <- read_csv('manuscript/brain/manuscript70/results/tables/subordinate_distruption_MEA_genes.csv')
head(sub)

sd <- cs[cs$symbol %in% sub$symbol,]
head(sd)

sd <- sd[!duplicated(sd[ , c("symbol")]),]


sd_upreg <- sd %>% filter(logFC <= -0.2)

sd_downreg <- sd %>% filter(logFC >= 0.2)


sd
# # we want the log2 fold change 
original_gene_list_up <- sd_upreg$logFC

original_gene_list_down <- sd_downreg$logFC

original_gene_list <- sd$logFC
# name the vector
names(original_gene_list_up) <- sd_upreg$symbol

names(original_gene_list_down) <- sd_downreg$symbol


names(original_gene_list) <- sd$symbol

# omit any NA values
up_list<-na.omit(original_gene_list_up)

down_list<-na.omit(original_gene_list_down)

gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
up_list = sort(up_list, decreasing = TRUE)

down_list = sort(down_list, decreasing = TRUE)

gene_list = sort(gene_list, decreasing = TRUE)



gse_up <- gseGO(geneList= up_list, 
                ont ="BP",
                keyType = 'SYMBOL', 
                nPerm = 10000,
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE, 
                OrgDb = organism, 
                pAdjustMethod = "none")

gse_down <- gseGO(geneList= down_list, 
                  ont ="BP",
                  keyType = 'SYMBOL', 
                  nPerm = 10000, 
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = organism, 
                  pAdjustMethod = "none")

gse <- gseGO(geneList= gene_list, 
             ont ="BP",
             keyType = 'SYMBOL', 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")




heatplot(gse_up, foldChange=up_list*-1, showCategory=5)
heatplot(gse_down, foldChange=down_list*-1, showCategory=5)
heatplot(gse, foldChange=gene_list*-1, showCategory=5)



