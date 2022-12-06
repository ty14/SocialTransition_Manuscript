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

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_mPFC_ControlDD.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol,ensgene, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

cdom <- limma_list$controldom
cdom$reg <- ifelse(cdom$logFC >= 0.2, "up", "down")

cdes <- limma_list$controldes 
cdes$reg <- ifelse(cdes$logFC >= 0.2, "up", "down")


domdes <- limma_list$domdes 
domdes$reg <- ifelse(domdes$logFC >= 0.2, "up", "down")



#looking for dom upregulated genes: 
cdd <- cdom %>% filter(reg == "down")

ddd <- domdes %>% filter(reg == "up")

up_dom <- cdd$symbol[(cdd$symbol %in% ddd$symbol)] %>% as.data.frame(.)
colnames(up_dom)<-"symbol"

up_dom<- up_dom %>% mutate(reg = "Up")

cddd <- cdom %>% filter(reg == "up")

dddd <- domdes %>% filter(reg == "down")

down_dom <- cddd$symbol[(cddd$symbol %in% dddd$symbol)] %>% as.data.frame(.)
colnames(down_dom)<-"symbol"

down_dom <- down_dom %>% mutate(reg = "Down")

#up dom = 52, down dom = 40
dom <- up_dom %>% full_join(down_dom)


write.csv(dom, 'manuscript/brain/manuscript70/results/tables/dominant_distruption_mPFC_genes.csv' ,row.names = F)


#looking further into descending genes:
dom <- read_csv('manuscript/brain/manuscript70/results/tables/dominant_distruption_mPFC_genes.csv')
head(dom)

cdom_top <- cdom[cdom$symbol %in% dom$symbol,]

cd_down <- cdom_top %>% arrange(-logFC)
x <- head(cd_down, 20)

dd_top <- domdes[domdes$symbol %in% dom$symbol,]
dd_down <- dd_top %>% arrange(logFC)
y <-head(dd_down, 20)
x[x$symbol %in% y$symbol,]


# UP: Krt12, Hkdc1, Ovol2, Grp,Rspo2, Ephb6,

#DOWN: Slc6a5, Gng10, Disc10, Ogg1,Dhrs13,Spata13 

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
    arrange(desc(GeneRatio)) %>% 
    mutate(direction = "Up") -> temp1
  
  
  fortify(
    ggo_down,
    showCategory = my_showCategory,
    by = "GeneRatio",
    split = NULL,
    includeAll = TRUE) %>% 
    arrange(desc(GeneRatio)) %>% 
    mutate(direction = "Down") -> temp2
  return(rbind(temp1,temp2))
  
}


gettop10GO(dom, my_showCategory) -> top10go1

write.csv(top10go1, 'manuscript/brain/manuscript70/results/tables/Dominant_mPFC_GOterms_BP.csv' ,row.names = F)


top10go1$Description[1:5]
top10go1$Description[6:10]
top10go1$Description[11:15]
top10go1$Description[16:20]
top10go1$Description[21:25]

up <- top10go1 %>% arrange(qvalue) %>% filter(direction == "Up")
head(up$Description)


down <- top10go1 %>% arrange(qvalue) %>% filter(direction == "Down")
head(down$Description)


########## cluster go terms stuff 

df <- read_csv('manuscript/brain/manuscript70/results/tables/Dominant_mPFC_GOterms_BP.csv')

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
df_upx.filename <-file.path("manuscript/brain/manuscript70/results/tables/mPFC_DDOM_up30_clusterprofiler_cluster_enr_results.txt")
write.table(df_upx,df_upx.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)



colnames(df_downx ) <- c("Name","Description", "pvalue","qvalue", "phenotype", "genes")
head(df_downx)
# write file for cytoscape
df_downx.filename <-file.path("manuscript/brain/manuscript70/results/tables/mPFC_DDOM_down30_clusterprofiler_cluster_enr_results.txt")
write.table(df_downx,df_downx.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)



### Make DGE boxplots:

ex <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_mPFC_ControlDD.RDS")

sm <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

sub_m_id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

# UP: Krt12, Hkdc1, Ovol2, Grp,Rspo2, Ephb6,


a1<- sm %>% filter(symbol == "Krt12")
a2 <- sm %>% filter(symbol == "Hkdc1")
a3 <- sm %>% filter(symbol == "Ovol2")
a4 <- sm %>% filter(symbol == "Grp")
a5<- sm %>% filter(symbol == "Rspo2")
a6 <- sm %>% filter(symbol == "Ephb6")


dp <- a1 %>% rbind(a2,a3,a4,a5,a6)


colnames(dp)

dp2 <- dp  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(sub_m_id)

dp2$group <- factor(dp2$group, levels = c("CDOM", "DOM", "DES")) 


source('functions/geom_boxjitter.R')


dp.l <- split(dp2 , dp2$symbol)
# UP: Krt12, *, *
# Hkdc1, *, **
# Ovol2, *, *
# Grp,*, *
# Rspo2,*,**
# Ephb6, *, *

c1 <- dp.l[[1]]
c1p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  scale_fill_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("mPFC Normalized Expression") +
  xlab("")+
  ggtitle(substitute(paste(italic("Ephb6"))))+
  ylim ( 3,9)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(8, 7),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))

c1p



c1 <- dp.l[[2]]
c2p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  scale_fill_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Grp"))))+
  ylim ( 3,8.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(7.5, 6.5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))


c2p



c1 <- dp.l[[3]]
c3p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  scale_fill_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Hkdc1"))))+
  ylim ( -2,7)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(6, 5.1),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "**"))

c3p


c1 <- dp.l[[6]]
c4p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  scale_fill_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Rspo2"))))+
  ylim ( 0,7.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(6.5, 5.5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "**"))

c4p


c1 <- dp.l[[5]]
c5p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  scale_fill_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Ovol2"))))+
  ylim ( 0,7.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(6.5, 5.5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))

c5p


c1 <- dp.l[[4]]
c6p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  scale_fill_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Krt12"))))+
  ylim ( -3,8)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(7,6),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))

c6p


domup <- gridExtra::grid.arrange(c1p, c2p, c5p, c4p, c3p, c6p, nrow = 1)

ggsave("manuscript/brain/manuscript70/results/results_figures/DOM_mPFC_TOPRegulated_genes_manuscriptSignifcants.png",domup, height = 3, width= 18,  dpi=300)

#DOWN: Slc6a5, Gng10, Disc10, Ogg1,Dhrs13,Spata13 


a1<- sm %>% filter(symbol == "Slc6a5")
a2 <- sm %>% filter(symbol == "Gng10")
a3 <- sm %>% filter(symbol == "Disc1")
a4 <- sm %>% filter(symbol == "Ogg1")
a5<- sm %>% filter(symbol == "Dhrs13")
a6 <- sm %>% filter(symbol == "Spata13")


dp <- a1 %>% rbind(a2,a3,a4,a5,a6)


colnames(dp)

dp2 <- dp  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(sub_m_id)

dp2$group <- factor(dp2$group, levels = c("CDOM", "DOM", "DES")) 


source('functions/geom_boxjitter.R')


dp.l <- split(dp2 , dp2$symbol)

#DOWN: Slc6a5, *, *
# Gng10, *, *
# Disc1, ***, **
# Ogg1, *, *
# Dhrs13, *, *
# Spata13 *, *


c1 <- dp.l[[1]]
c1p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  scale_fill_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab(" MEA Normalized Expression") +
  xlab("")+
  ggtitle(substitute(paste(italic("Dhrs13"))))+
  ylim ( 0 ,6)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(5, 4),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))


c1p



c1 <- dp.l[[2]]
c2p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  scale_fill_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Disc1"))))+
  ylim ( 1,7)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(6, 5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("***", "**"))


c2p

c1 <- dp.l[[3]]
c3p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  scale_fill_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Gng10"))))+
  ylim ( -3,6)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(5, 4),
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
  scale_color_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  scale_fill_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Ogg1"))))+
  ylim ( 0, 5.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(4.5,3.5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))


c4p



c1 <- dp.l[[5]]
c5p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  scale_fill_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Slc6a5"))))+
  ylim ( -1,6.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(5.5, 4.5),
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
  scale_color_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  scale_fill_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Spata13"))))+
  ylim ( 1,7)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(6, 5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))


c6p



domdown <- gridExtra::grid.arrange(c1p, c2p, c3p, c4p, c5p, c6p, nrow = 1)




ggsave("manuscript/brain/manuscript70/results/results_figures/DOM_mPFC_DOWNRegulated_genes_manuscriptSignifcants.png",domdown, height = 3, width= 18,  dpi=300)


##################################################
#Ascenders

my_logFC_threshold = 0.2

#mPFC data

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_mPFC_ControlSUB.RDS") %>% 
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

#up dom = 89, down dom = 103
sub <- up_sub%>% full_join(down_sub)



write.csv(sub, 'manuscript/brain/manuscript70/results/tables/subordinate_distruption_mPFC_genes.csv' ,row.names = F)

#looking further into staying subordinate  genes:
sub<- read_csv('manuscript/brain/manuscript70/results/tables/subordinate_distruption_mPFC_genes.csv')
head(sub)


casc_top <- cs[cs$symbol %in% sub$symbol,]

cd_up <- casc_top %>% arrange(-logFC)
x <- head(cd_up, 15)

sa_top <- sa[sa$symbol %in% sub$symbol,]
sa_up <- sa_top %>% arrange(logFC)
y <- head(sa_up, 15)


x[x$symbol %in% y$symbol,]

#UP: St8sia2, Vmn2r85, Etnk2, Hrh2, Iigp1, Cracr2a

#DOWN:  Arhgap36, Magel2, Ndst4, Plekha2, Dlk1, Cntnap5c 




gettop10GO(sub, my_showCategory) -> top10go1

write.csv(top10go1, 'manuscript/brain/manuscript70/results/tables/SUB_mPFC_GOterms_BP.csv' ,row.names = F)


top10go1$Description[1:5]
top10go1$Description[6:10]
top10go1$Description[11:15]
top10go1$Description[16:20]
top10go1$Description[21:25]

up <- top10go1 %>% arrange(qvalue) %>% filter(direction == "Up")
head(up$Description)


down <- top10go1 %>% arrange(qvalue) %>% filter(direction == "Down")
head(down$Description)


########## cluster go terms stuff 

df <- read_csv('manuscript/brain/manuscript70/results/tables/SUB_mPFC_GOterms_BP.csv')

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
df_upx.filename <-file.path("manuscript/brain/manuscript70/results/tables/mPFC_SUB_up30_clusterprofiler_cluster_enr_results.txt")
write.table(df_upx,df_upx.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)



colnames(df_downx ) <- c("Name","Description", "pvalue","qvalue", "phenotype", "genes")
head(df_downx)
# write file for cytoscape
df_downx.filename <-file.path("manuscript/brain/manuscript70/results/tables/mPFC_SUB_down30_clusterprofiler_cluster_enr_results.txt")
write.table(df_downx,df_downx.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)


### Make DGE boxplots:

ex <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_mPFC_ControlSUB.RDS")

sm <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

sub_m_id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

# St8sia2, Vmn2r85, Etnk2, Hrh2, Iigp1, Cracr2a

a1<- sm %>% filter(symbol == "St8sia2")
a2 <- sm %>% filter(symbol == "Vmn2r85")
a3 <- sm %>% filter(symbol == "Etnk2")
a4 <- sm %>% filter(symbol == "Hrh2")
a5<- sm %>% filter(symbol == "Iigp1")
a6 <- sm %>% filter(symbol == "Cracr2a")


dp <- a1 %>% rbind(a2,a3,a4,a5,a6)


colnames(dp)

dp2 <- dp  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(sub_m_id)

dp2$group <- factor(dp2$group, levels = c("CSUB", "SUB", "ASC")) 


source('functions/geom_boxjitter.R')


dp.l <- split(dp2 , dp2$symbol)
# St8sia2, ***, *
# Vmn2r85, **, * 
# Etnk2, *, **
# Hrh2, *, *
# Iigp1,*, *
# Cracr2a **, **
# MEA Normalized Expression

c1 <- dp.l[[5]]
c1p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  scale_fill_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("mPFC Normalized Expression") +
  xlab("")+
  ggtitle(substitute(paste(italic("St8sia2"))))+
  ylim (0,7)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(6, 5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("***", "*"))

c1p


c1 <- dp.l[[3]]
c2p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  scale_fill_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Hrh2"))))+
  ylim ( 2,7.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(7, 6),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))

c2p





c1 <- dp.l[[1]]
c3p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  scale_fill_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Cracr2a"))))+
  ylim ( 1,6.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(5.5,4.5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "**"))

c3p


c1 <- dp.l[[4]]
c4p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  scale_fill_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Iigp1"))))+
  ylim ( 0,5)+
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

c4p


c1 <- dp.l[[2]]
c5p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  scale_fill_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Etnk2"))))+
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
              annotations = c("*", "**"))

c5p


c1 <- dp.l[[6]]
c6p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  scale_fill_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Vmn2r85"))))+
  ylim ( -1, 5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(4, 3),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "*"))

c6p


subup <- gridExtra::grid.arrange(c1p, c2p, c5p, c4p, c3p, c6p, nrow = 1)

ggsave("manuscript/brain/manuscript70/results/results_figures/SUB_mPFC_TOPRegulated_genes_manuscriptSignifcants.png",subup, height = 3, width= 18,  dpi=300)

#  Arhgap36, Magel, Ndst4, Plekha2, Dlk1, Cntnap5c 


a1<- sm %>% filter(symbol == "Arhgap36")
a2 <- sm %>% filter(symbol == "Magel2")
a3 <- sm %>% filter(symbol == "Ndst4")
a4 <- sm %>% filter(symbol == "Plekha2")
a5<- sm %>% filter(symbol == "Dlk1")
a6 <- sm %>% filter(symbol == "Cntnap5c")


dp <- a1 %>% rbind(a2,a3,a4,a5,a6)


colnames(dp)

dp2 <- dp  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(sub_m_id)

dp2$group <- factor(dp2$group, levels = c("CSUB", "SUB", "ASC")) 


source('functions/geom_boxjitter.R')


dp.l <- split(dp2 , dp2$symbol)
# Arhgap36, **, *
# Magel, **, *
# Ndst4,*, *
# Plekha2, ***, **
# Dlk1, *, *
# Cntnap5c *, *


c1 <- dp.l[[5]]
c1p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  scale_fill_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("mPFC Normalized Expression") +
  xlab("")+
  ggtitle(substitute(paste(italic("Ndst4"))))+
  ylim (1,8.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(7.5, 6.5),
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
  scale_color_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  scale_fill_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Dlk1"))))+
  ylim ( 0,7)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(6, 5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))

c2p





c1 <- dp.l[[1]]
c3p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  scale_fill_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Arhgap36"))))+
  ylim ( -3,5.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(4.5, 3.5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "*"))

c3p


c1 <- dp.l[[4]]
c4p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  scale_fill_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Magel2"))))+
  ylim ( -3,5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(3.75, 2.75),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "*"))

c4p


c1 <- dp.l[[2]]
c5p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  scale_fill_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Cntnap5c"))))+
  ylim ( -3, 4)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(3, 2),
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
  scale_color_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  scale_fill_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Plekha2"))))+
  ylim ( 2, 8)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "SUB"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(7, 6),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("***", "**"))

c6p


subdown <- gridExtra::grid.arrange(c1p, c2p, c5p, c4p, c3p, c6p, nrow = 1)

ggsave("manuscript/brain/manuscript70/results/results_figures/SUB_DOWN_mPFC_TOPRegulated_genes_manuscriptSignifcants.png",subdown, height = 3, width= 18,  dpi=300)


