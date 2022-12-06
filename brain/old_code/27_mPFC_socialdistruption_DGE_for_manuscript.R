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



#looking for cdom upregulated genes: 
cdd <- cdom %>% filter(reg == "up")

ddd <- cdes %>% filter(reg == "up")

up_dom <- cdd$symbol[(cdd$symbol %in% ddd$symbol)] %>% as.data.frame(.)
colnames(up_dom)<-"symbol"

up_dom<- up_dom %>% mutate(reg = "Up")

cddd <- cdom %>% filter(reg == "down")

dddd <- cdes %>% filter(reg == "down")

down_dom <- cddd$symbol[(cddd$symbol %in% dddd$symbol)] %>% as.data.frame(.)
colnames(down_dom)<-"symbol"

down_dom <- down_dom %>% mutate(reg = "Down")

#up dom = 96, down dom = 119
dom <- up_dom %>% full_join(down_dom)

write.csv(dom, 'manuscript/brain/manuscript70/results/tables/control_dominant_distruption_mPFC_genes.csv' ,row.names = F)


#looking further into descending genes:
dom <- read_csv('manuscript/brain/manuscript70/results/tables/control_dominant_distruption_mPFC_genes.csv')
head(dom)



cdom_top <- cdom[cdom$symbol %in% dom$symbol,]

cd_down <- cdom_top %>% arrange(logFC)
x <- head(cd_down, 20)

dd_top <- cdes[cdes$symbol %in% dom$symbol,]
dd_down <- dd_top %>% arrange(logFC)
y <-head(dd_down, 20)
x[x$symbol %in% y$symbol,]


# UP: Dmkn, Ankrd66, Gm867, Aqp1, Eomes, Gbx2

#DOWN: Cd14, Csrnp1, Egr4, Fosb, Calhm2, Sparr



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

write.csv(top10go1, 'manuscript/brain/manuscript70/results/tables/CDOM_mPFC_GOterms_BP.csv' ,row.names = F)


########## cluster go terms stuff 

df <- read_csv('manuscript/brain/manuscript70/results/tables/CDOM_mPFC_GOterms_BP.csv')

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
df_upx.filename <-file.path("manuscript/brain/manuscript70/results/tables/mPFC_CDOM_up30_clusterprofiler_cluster_enr_results.txt")
write.table(df_upx,df_upx.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)



colnames(df_downx ) <- c("Name","Description", "pvalue","qvalue", "phenotype", "genes")
head(df_downx)
# write file for cytoscape
df_downx.filename <-file.path("manuscript/brain/manuscript70/results/tables/mPFC_CDOM_down30_clusterprofiler_cluster_enr_results.txt")
write.table(df_downx,df_downx.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)




### Make DGE boxplots:

ex <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_mPFC_ControlDD.RDS")

sm <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

sub_m_id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

# UP: Dmkn, Ankrd66, Gm867, Aqp1, Eomes, Gbx2


a1<- sm %>% filter(symbol == "Dmkn")
a2 <- sm %>% filter(symbol == "Ankrd66")
a3 <- sm %>% filter(symbol == "Gm867")
a4 <- sm %>% filter(symbol == "Aqp1")
a5<- sm %>% filter(symbol == "Eomes")
a6 <- sm %>% filter(symbol == "Gbx2")


dp <- a1 %>% rbind(a2,a3,a4,a5,a6)


colnames(dp)

dp2 <- dp  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(sub_m_id)

dp2$group <- factor(dp2$group, levels = c("CDOM", "DOM", "DES")) 


source('functions/geom_boxjitter.R')


dp.l <- split(dp2 , dp2$symbol)

# Dmkn, **, **
# Ankrd66, *, *
# Gm867, **, **
# Aqp1, *, *
# Eomes, *,*
# Gbx2 *, *


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
  ggtitle(substitute(paste(italic("Ankrd66"))))+
  ylim ( -3,4)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("CDOM", "DOM")),
              map_signif_level = TRUE,
              y_position = c(3.4, 2.5),
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
  ggtitle(substitute(paste(italic("Aqp1"))))+
  ylim ( -3,4.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("CDOM", "DOM")),
              map_signif_level = TRUE,
              y_position = c(3.75, 2.90),
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
  ylim( 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Dmkn"))))+
  ylim ( -3.5, 5.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("CDOM", "DOM")),
              map_signif_level = TRUE,
              y_position = c(4.5, 3.5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "**"))

c3p


c1 <- dp.l[[5]]
c4p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  scale_fill_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  ylim( 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Gbx2"))))+
  ylim ( -3,5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("CDOM", "DOM")),
              map_signif_level = TRUE,
              y_position = c(4, 3),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))

c4p


c1 <- dp.l[[6]]
c5p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  scale_fill_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Gm867"))))+
  ylim ( -3,5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("CDOM", "DOM")),
              map_signif_level = TRUE,
              y_position = c(4.25, 3.25),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "**"))

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
  ggtitle(substitute(paste(italic("Eomes"))))+
  ylim ( -1,8.25)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("CDOM", "DOM")),
              map_signif_level = TRUE,
              y_position = c(7.25,6.25),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))

c6p


domup <- gridExtra::grid.arrange(c1p, c2p, c5p, c4p, c3p, c6p, nrow = 1)

ggsave("manuscript/brain/manuscript70/results/results_figures/CDOM_mPFC_TOPRegulated_genes_manuscriptSignifcants.png",domup, height = 3, width= 18,  dpi=300)


# Cd14, Csrnp1, Egr4, Fosb, Calhm2, Spaar


a1<- sm %>% filter(symbol == "Cd14")
a2 <- sm %>% filter(symbol == "Csrnp1")
a3 <- sm %>% filter(symbol == "Egr4")
a4 <- sm %>% filter(symbol == "Fosb")
a5<- sm %>% filter(symbol == "Calhm2")
a6 <- sm %>% filter(symbol == "Spaar")


dp <- a1 %>% rbind(a2,a3,a4,a5,a6)


colnames(dp)

dp2 <- dp  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(sub_m_id)

dp2$group <- factor(dp2$group, levels = c("CDOM", "DOM", "DES")) 


source('functions/geom_boxjitter.R')


dp.l <- split(dp2 , dp2$symbol)

# Cd14, **, *
# Csrnp1,**, ***
# Egr4, **, *
# Fosb, *, *
# Calhm2, *, **
# Spaar *, *


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
  ggtitle(substitute(paste(italic("Calhm2"))))+
  ylim ( 0,6)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("CDOM", "DOM")),
              map_signif_level = TRUE,
              y_position = c(5, 4),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "**"))

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
  ggtitle(substitute(paste(italic("Cd14"))))+
  ylim ( -3,4.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("CDOM", "DOM")),
              map_signif_level = TRUE,
              y_position = c(3.75, 2.90),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "*"))

c2p



c1 <- dp.l[[3]]
c3p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  scale_fill_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  ylim( 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Csrnp1"))))+
  ylim ( -1, 6.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("CDOM", "DOM")),
              map_signif_level = TRUE,
              y_position = c(5.5, 4.5),
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
  ylim( 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Spaar"))))+
  ylim ( -1,5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("CDOM", "DOM")),
              map_signif_level = TRUE,
              y_position = c(4, 3),
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
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Fosb"))))+
  ylim ( -1,5.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("CDOM", "DOM")),
              map_signif_level = TRUE,
              y_position = c(4.5, 3.5),
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
  ggtitle(substitute(paste(italic("Egr4"))))+
  ylim ( 2,10)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("CDOM", "DOM")),
              map_signif_level = TRUE,
              y_position = c(9,8),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "*"))

c6p


domdown <- gridExtra::grid.arrange(c1p, c2p, c5p, c4p, c3p, c6p, nrow = 1)

ggsave("manuscript/brain/manuscript70/results/results_figures/CDOM_mPFC_DOWN_TOPRegulated_genes_manuscriptSignifcants.png",domdown, height = 3, width= 18,  dpi=300)




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

#looking for csub upregulated genes: 
cdd <- cs %>% filter(reg == "up")

ddd <- ca %>% filter(reg == "up")

up_sub <- cdd$symbol[(cdd$symbol %in% ddd$symbol)] %>% as.data.frame(.)
colnames(up_sub)<-"symbol"

up_sub<- up_sub %>% mutate(reg = "Up")

cddd <- cs %>% filter(reg == "down")

dddd <- ca %>% filter(reg == "down")

down_sub <- cddd$symbol[(cddd$symbol %in% dddd$symbol)] %>% as.data.frame(.)
colnames(down_sub)<-"symbol"

down_sub <- down_sub %>% mutate(reg = "Down")

#up sub = 120, down sub = 98
sub <- up_sub %>% full_join(down_sub)

write.csv(sub, 'manuscript/brain/manuscript70/results/tables/control_subordinate_distruption_mPFC_genes.csv' ,row.names = F)


#looking further into descending genes:
sub <- read_csv('manuscript/brain/manuscript70/results/tables/control_subordinate_distruption_mPFC_genes.csv')
head(sub)



csub_top <- cs[cs$symbol %in% sub$symbol,]

cd_down <- csub_top %>% arrange(-logFC)
x <- head(cd_down, 20)

dd_top <- ca[ca$symbol %in% sub$symbol,]
dd_down <- dd_top %>% arrange(-logFC)
y <-head(dd_down, 20)
x[x$symbol %in% y$symbol,]


# UP: Vmo1, Serpina3i, Ccl12, Vstm4, Prdm12, Crabp1

#DOWN: Tmem232, Caps2, Scn7a, Cbln1, Opn3, Zfp804b



gettop10GO(sub, my_showCategory) -> top10go1

write.csv(top10go1, 'manuscript/brain/manuscript70/results/tables/CSUB_mPFC_GOterms_BP.csv' ,row.names = F)


########## cluster go terms stuff 

df <- read_csv('manuscript/brain/manuscript70/results/tables/CSUB_mPFC_GOterms_BP.csv')

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
df_upx.filename <-file.path("manuscript/brain/manuscript70/results/tables/mPFC_CSUB_up30_clusterprofiler_cluster_enr_results.txt")
write.table(df_upx,df_upx.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)



colnames(df_downx ) <- c("Name","Description", "pvalue","qvalue", "phenotype", "genes")
head(df_downx)
# write file for cytoscape
df_downx.filename <-file.path("manuscript/brain/manuscript70/results/tables/mPFC_CSUB_down30_clusterprofiler_cluster_enr_results.txt")
write.table(df_downx,df_downx.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)



### Make DGE boxplots:

ex <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_mPFC_ControlSUB.RDS")

sm <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

sub_m_id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

# Vmo1, Serpina3i, Ccl12, Vstm4, Prdm12, Crabp1

a1<- sm %>% filter(symbol == "Serpina3i")
a2 <- sm %>% filter(symbol == "Vmo1")
a3 <- sm %>% filter(symbol == "Ccl12")
a4 <- sm %>% filter(symbol == "Vstm4")
a5<- sm %>% filter(symbol == "Prdm12")
a6 <- sm %>% filter(symbol == "Crabp1")


dp <- a1 %>% rbind(a2,a3,a4,a5,a6)


colnames(dp)

dp2 <- dp  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(sub_m_id)

dp2$group <- factor(dp2$group, levels = c("CSUB", "SUB", "ASC")) 


source('functions/geom_boxjitter.R')


dp.l <- split(dp2 , dp2$symbol)
# Vmo1, *, ***
# Serpina3i, *,*
# Ccl12, **, **
# Vstm4, *, *
# Prdm12, *,*
# Crabp1 *, *

c1 <- dp.l[[4]]
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
  ggtitle(substitute(paste(italic("Serpina3i"))))+
  ylim (-3,6)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("CSUB", "SUB")),
              map_signif_level = TRUE,
              y_position = c(5, 4),
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
  ggtitle(substitute(paste(italic("Prdm12"))))+
  ylim (-3,5.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("CSUB", "SUB")),
              map_signif_level = TRUE,
              y_position = c(4.75, 3.75),
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
  ggtitle(substitute(paste(italic("Ccl12"))))+
  ylim (-3,5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("CSUB", "SUB")),
              map_signif_level = TRUE,
              y_position = c(4,3),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "**"))

c3p


c1 <- dp.l[[5]]
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
  ggtitle(substitute(paste(italic("Vmo1"))))+
  ylim ( -3,5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("CSUB", "SUB")),
              map_signif_level = TRUE,
              y_position = c(4, 3),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "***"))

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
  ggtitle(substitute(paste(italic("Crabp1"))))+
  ylim ( -3,5.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("CSUB", "SUB")),
              map_signif_level = TRUE,
              y_position = c(4.5, 3.5),
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
  ggtitle(substitute(paste(italic("Vstm4"))))+
  ylim ( -3, 6)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("CSUB", "SUB")),
              map_signif_level = TRUE,
              y_position = c(5, 4),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))

c6p


subup <- gridExtra::grid.arrange(c1p, c2p, c5p, c4p, c3p, c6p, nrow = 1)

ggsave("manuscript/brain/manuscript70/results/results_figures/CSUB_mPFC_TOPRegulated_genes_manuscriptSignifcants.png",domup, height = 3, width= 18,  dpi=300)


# Tmem232, Caps2, Scn7a, Cbln1, Opn3, Zfp804b

a1<- sm %>% filter(symbol == "Tmem232")
a2 <- sm %>% filter(symbol == "Caps2")
a3 <- sm %>% filter(symbol == "Scn7a")
a4 <- sm %>% filter(symbol == "Opn3")
a5<- sm %>% filter(symbol == "Cbln1")
a6 <- sm %>% filter(symbol == "Zfp804b")


dp <- a1 %>% rbind(a2,a3,a4,a5,a6)


colnames(dp)

dp2 <- dp  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(sub_m_id)

dp2$group <- factor(dp2$group, levels = c("CSUB", "SUB", "ASC")) 


source('functions/geom_boxjitter.R')


dp.l <- split(dp2 , dp2$symbol)
# Tmem232, **,***
# Caps2, *, **
# Scn7a, *, **
# Cbln1, **, *
# Opn3, *, **
# Zfp804b *,**

c1 <- dp.l[[4]]
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
  ggtitle(substitute(paste(italic("Scn7a"))))+
  ylim (0,7)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("CSUB", "SUB")),
              map_signif_level = TRUE,
              y_position = c(6, 5),
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
  scale_color_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  scale_fill_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Opn3"))))+
  ylim (-0,6)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("CSUB", "SUB")),
              map_signif_level = TRUE,
              y_position = c(5.25,4.25),
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
  scale_color_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  scale_fill_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Caps2"))))+
  ylim (-1,6)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("CSUB", "SUB")),
              map_signif_level = TRUE,
              y_position = c(5,4),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "**"))

c3p


c1 <- dp.l[[5]]
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
  ggtitle(substitute(paste(italic("Tmem232"))))+
  ylim ( -3,5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("CSUB", "SUB")),
              map_signif_level = TRUE,
              y_position = c(4, 3),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "***"))

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
  ggtitle(substitute(paste(italic("Cbln1"))))+
  ylim ( -2,5.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("CSUB", "SUB")),
              map_signif_level = TRUE,
              y_position = c(4.5, 3.5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "*"))

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
  ggtitle(substitute(paste(italic("Zfp804b"))))+
  ylim ( -1, 6)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("CSUB", "SUB")),
              map_signif_level = TRUE,
              y_position = c(5, 4),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "**"))

c6p


subup <- gridExtra::grid.arrange(c1p, c2p, c5p, c4p, c3p, c6p, nrow = 1)

ggsave("manuscript/brain/manuscript70/results/results_figures/CSUB_mPFC_DOWN_TOPRegulated_genes_manuscriptSignifcants.png",subup, height = 3, width= 18,  dpi=300)
