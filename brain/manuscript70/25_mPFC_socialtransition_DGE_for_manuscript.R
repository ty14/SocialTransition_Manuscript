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

#mPFC data

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
cdd <- cdes %>% filter(reg == "down")

ddd <- domdes %>% filter(reg == "down")

up_dom <- cdd$symbol[(cdd$symbol %in% ddd$symbol)] %>% as.data.frame(.)
colnames(up_dom)<-"symbol"

up_dom<- up_dom %>% mutate(reg = "Up")

cddd <- cdes %>% filter(reg == "up")

dddd <- domdes %>% filter(reg == "up")

down_dom <- cddd$symbol[(cddd$symbol %in% dddd$symbol)] %>% as.data.frame(.)
colnames(down_dom)<-"symbol"

down_dom <- down_dom %>% mutate(reg = "Down")

#up sd= 85 , down dom = 62
des <- up_dom %>% full_join(down_dom)

write.csv(des, 'manuscript/brain/manuscript70/results/tables/descenders_tranisiton_mPFC_genes.csv' ,row.names = F)


#looking further into descending genes:
dom <- read_csv('manuscript/brain/manuscript70/results/tables/descenders_tranisiton_mPFC_genes.csv')
head(dom)



cdes_top <- cdes[cdes$symbol %in% dom$symbol,]

cd_down <- cdes_top %>% arrange(-logFC)
x <- head(cd_down, 20)

dd_top <- domdes[domdes$symbol %in% dom$symbol,]
dd_down <- dd_top %>% arrange(-logFC)
y <-head(dd_down, 20)
x[x$symbol %in% y$symbol,]


# UP: S100a9, Ngp, Wfdc18, Ccl6, Olfml2b, Htra4

#DOWN: Nkpd1, Gchfr,Erbb3, Morc1, Ntng1, Mak



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

write.csv(top10go1, 'manuscript/brain/manuscript70/results/tables/SD_mPFC_GOterms_BP.csv' ,row.names = F)

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

df <- read_csv('manuscript/brain/manuscript70/results/tables/SD_mPFC_GOterms_BP.csv')

head(df) 

#up data 
df_up <- df %>%  filter(direction == "Up") 
df_upx <- df_up[1:21,]

# down data 
df_down <- df %>%  filter(direction == "Down") 
df_downx <- df_down[1:21,]

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
df_upx.filename <-file.path("manuscript/brain/manuscript70/results/tables/mPFC_SD_up30_clusterprofiler_cluster_enr_results.txt")
write.table(df_upx,df_upx.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)



colnames(df_downx ) <- c("Name","Description", "pvalue","qvalue", "phenotype", "genes")
head(df_downx)
# write file for cytoscape
df_downx.filename <-file.path("manuscript/brain/manuscript70/results/tables/mPFC_SD_down30_clusterprofiler_cluster_enr_results.txt")
write.table(df_downx,df_downx.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)





### Make DGE boxplots:

ex <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_mPFC_ControlDD.RDS")

sm <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

sub_m_id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

# UP: S100a9, Ngp, Wfdc18, Ccl6, Olfml2b, Htra4

t <- sm %>% filter(symbol == "Thg1l")

a1<- sm %>% filter(symbol == "S100a9")
a2 <- sm %>% filter(symbol == "Ngp")
a3 <- sm %>% filter(symbol == "Wfdc18")
a4 <- sm %>% filter(symbol == "Ccl6")
a5<- sm %>% filter(symbol == "Olfml2b")
a6 <- sm %>% filter(symbol == "Htra4")

ma <- sm %>% filter(symbol == "Maoa")
mb <- sm %>% filter(symbol == "Maob")


dp <- a1 %>% rbind(a2,a3,a4,a5,a6)

dp <- ma %>% rbind(mb)
colnames(dp)

dp2 <- t  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(sub_m_id)

dp2$group <- factor(dp2$group, levels = c("CDOM", "DOM", "DES")) 


source('functions/geom_boxjitter.R')


dp.l <- split(dp2 , dp2$symbol)

# S100a9, *, *
# Ngp, *, *
# Wfdc18, **, **
# Ccl6, *,*
# Olfml2b, **, **
# Htra4 *, **
# MEA Normalized Expression


c1 <- dp.l[[2]]
c1p <-  ggplot(dp2, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  scale_fill_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("mPFC Normalized Expression") +
  xlab("")+
  ggtitle(substitute(paste(italic("Thg1l"))))+
  ylim (1, 6)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(5.25, 4.25),
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
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Htra4"))))+
  ylim ( -3,5.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(4.5, 3.5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "**"))

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
  ggtitle(substitute(paste(italic("Ngp"))))+
  ylim ( -3,5.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(4.5, 3.5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))

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
  ggtitle(substitute(paste(italic("Wfdc18"))))+
  ylim ( 0,6.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(5.75, 4.75),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "**"))

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
  ggtitle(substitute(paste(italic("S100a9"))))+
  ylim ( -3,9.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(8.5, 7.25),
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
  ggtitle(substitute(paste(italic("Olfml2b"))))+
  ylim ( -2.5,6)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(5,4),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "**"))

c6p


domup <- gridExtra::grid.arrange(c1p, c2p, c5p, c4p, c3p, c6p, nrow = 1)


ggsave("manuscript/brain/manuscript70/results/results_figures/SD_mPFC_TOPRegulated_genes_manuscriptSignifcants.png",domup, height = 3, width= 18,  dpi=300)



## DOWN

a1<- sm %>% filter(symbol == "Nkpd1")
a2 <- sm %>% filter(symbol == "Gchfr")
a3 <- sm %>% filter(symbol == "Erbb3")
a4 <- sm %>% filter(symbol == "Morc1")
a5<- sm %>% filter(symbol == "Ntng1")
a6 <- sm %>% filter(symbol == "Mak")



table(dp2$symbol)
colnames(dp)

dp2 <- dp  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(sub_m_id)

dp2$group <- factor(dp2$group, levels = c("CDOM", "DOM", "DES")) 


source('functions/geom_boxjitter.R')


dp.l <- split(dp2 , dp2$symbol)




#DOWN:  Nkpd1, **, **
# Gchfr, *, **
# Erbb3, *, *
# Morc1, ***, *
# Ntng1, **, *
# Mak *, *

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
  ggtitle(substitute(paste(italic("Erbb3"))))+
  ylim ( -2,5.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(4.5, 3.5),
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
  ggtitle(substitute(paste(italic("Gchfr"))))+
  ylim ( -2.5, 5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(4, 3),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "**"))

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
  ggtitle(substitute(paste(italic("Mak"))))+
  ylim ( -0,5.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(4.5, 3.5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))

c3p


c1 <- dp.l[[5]]
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
  ggtitle(substitute(paste(italic("Nkpd1"))))+
  ylim ( -3,4.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(3.5, 2.5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "**"))

c4p


c1 <- dp.l[[6]]
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
  ggtitle(substitute(paste(italic("Ntng1"))))+
  ylim ( 3,8.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(7.5, 6.5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "*"))

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
  ggtitle(substitute(paste(italic("Morc1"))))+
  ylim ( 1.5,7)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DES"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(6,5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("***", "*"))

c6p


domdown <- gridExtra::grid.arrange(c1p, c2p, c3p, c4p, c5p, c6p, nrow = 1)




ggsave("manuscript/brain/manuscript70/results/results_figures/SD_mPFC_DOWNRegulated_genes_manuscriptSignifcants.png",domdown, height = 3, width= 18,  dpi=300)




###############################################
# Ascenders
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
cdd <- ca %>% filter(reg == "down")

ddd <- sa %>% filter(reg == "down")

up_sub <- cdd$symbol[(cdd$symbol %in% ddd$symbol)] %>% as.data.frame(.)
colnames(up_sub)<-"symbol"

up_sub<- up_sub %>% mutate(reg = "Up")

cddd <- ca %>% filter(reg == "up")

dddd <- sa %>% filter(reg == "up")

down_sub <- cddd$symbol[(cddd$symbol %in% dddd$symbol)] %>% as.data.frame(.)
colnames(down_sub)<-"symbol"

down_sub <- down_sub %>% mutate(reg = "Down")

#up dom = 58, down dom = 60
asc <- up_sub%>% full_join(down_sub)



write.csv(asc, 'manuscript/brain/manuscript70/results/tables/ASC_transition_mPFC_genes.csv' ,row.names = F)

#looking further into staying subordinate  genes:
asc<- read_csv('manuscript/brain/manuscript70/results/tables/ASC_transition_mPFC_genes.csv')
head(sub)


casc_top <- ca[ca$symbol %in% asc$symbol,]

cd_up <- casc_top %>% arrange(logFC)
x <- head(cd_up, 15)

sa_top <- sa[sa$symbol %in% asc$symbol,]
sa_up <- sa_top %>% arrange(logFC)
y <- head(sa_up, 15)

x[x$symbol %in% y$symbol,]

#UP: Ttr, Cfap100, Slc25a34, Gtpbp10, Eme2, Afap1l1

#DOWN:  Lrrc25, Sell, Fermt1, Mpeg1, Six6, Olfr127


gettop10GO(asc, my_showCategory) -> top10go1

write.csv(top10go1, 'manuscript/brain/manuscript70/results/tables/ASC_mPFC_GOterms_BP.csv' ,row.names = F)

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

df <- read_csv('manuscript/brain/manuscript70/results/tables/ASC_mPFC_GOterms_BP.csv')

head(df) 

#up data 
df_up <- df %>%  filter(direction == "Up") 
df_upx <- df_up[1:13,]

# down data 
df_down <- df %>%  filter(direction == "Down") 
df_downx <- df_down[1:13,]

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
df_upx.filename <-file.path("manuscript/brain/manuscript70/results/tables/mPFC_ASC_up30_clusterprofiler_cluster_enr_results.txt")
write.table(df_upx,df_upx.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)



colnames(df_downx ) <- c("Name","Description", "pvalue","qvalue", "phenotype", "genes")
head(df_downx)
# write file for cytoscape
df_downx.filename <-file.path("manuscript/brain/manuscript70/results/tables/mPFC_ASC_down30_clusterprofiler_cluster_enr_results.txt")
write.table(df_downx,df_downx.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)


### Make DGE boxplots:

ex <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_mPFC_ControlSUB.RDS")

sm <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

sub_m_id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

#  Ttr, Cfap100, Slc25a34, Gtpbp10, Eme2, Afap1l1

a1<- sm %>% filter(symbol == "Ttr")
a2 <- sm %>% filter(symbol == "Cfap100")
a3 <- sm %>% filter(symbol == "Slc25a34")
a4 <- sm %>% filter(symbol == "Gtpbp10")
a5<- sm %>% filter(symbol == "Eme2")
a6 <- sm %>% filter(symbol == "Afap1l1")

t <- sm %>% filter(symbol == "Thg1l")


dp <- a1 %>% rbind(a2,a3,a4,a5,a6)


ma <- sm %>% filter(symbol == "Maob")
colnames(dp)

dp2 <- t%>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(sub_m_id)

dp2$group <- factor(dp2$group, levels = c("CSUB", "SUB", "ASC")) 


source('functions/geom_boxjitter.R')


dp.l <- split(dp2 , dp2$symbol)
# Ttr, *, **
# Cfap100, **, *
# Slc25a34, *, *
# Gtpbp10, *, *
#Eme2, *, ***
# Afap1l1 *, **

c1 <- dp.l[[5]]
c1p <-  ggplot(dp2, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  scale_fill_manual(values = c("#f9b641ff", "#b8627dff", "#414487FF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("mPFC Normalized Expression") +
  xlab("")+
  ggtitle(substitute(paste(italic("Thg1l"))))+
  ylim (1.5,5.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("SUB", "ASC")),
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
  ggtitle(substitute(paste(italic("Eme2"))))+
  ylim ( 2,6)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(5.5, 4.5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "***"))

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
  ggtitle(substitute(paste(italic(" Afap1l1"))))+
  ylim ( 2,8)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(7,6),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "**"))

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
  ggtitle(substitute(paste(italic("Gtpbp10"))))+
  ylim ( -0,5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(4.25, 3.25),
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
  ggtitle(substitute(paste(italic("Cfap100"))))+
  ylim ( -0,6)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(5, 4),
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
  ggtitle(substitute(paste(italic("Ttr"))))+
  ylim ( 2.5, 10.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", 'ASC'),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(9.5, 8.5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "**"))

c6p


subup <- gridExtra::grid.arrange(c1p, c2p, c5p, c4p, c3p, c6p, nrow = 1)

ggsave("manuscript/brain/manuscript70/results/results_figures/ASC_mPFC_TOPRegulated_genes_manuscriptSignifcants.png",subup, height = 3, width= 18,  dpi=300)


#  DOWN: Lrrc25, Sell, Fermt1, Mpeg1, Six6, Olfr127
a1<- sm %>% filter(symbol == "Lrrc25")
a2 <- sm %>% filter(symbol == "Sell")
a3 <- sm %>% filter(symbol == "Fermt1")
a4 <- sm %>% filter(symbol == "Mpeg1")
a5<- sm %>% filter(symbol == "Six6")
a6 <- sm %>% filter(symbol == "Olfr127")


dp <- a1 %>% rbind(a2,a3,a4,a5,a6)


colnames(dp)

dp2 <- dp  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(sub_m_id)

dp2$group <- factor(dp2$group, levels = c("CSUB", "SUB", "ASC")) 


source('functions/geom_boxjitter.R')


dp.l <- split(dp2 , dp2$symbol)
# Lrrc25, **, *** 
# Sell, **, **
# Fermt1, **, *
# Mpeg1, *, *
# Six6, *, *
# Olfr127 **, ***

c1 <- dp.l[[1]]
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
  ggtitle(substitute(paste(italic("Fermt1"))))+
  ylim (-2,4)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(3.25, 2.25),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "*"))

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
  ggtitle(substitute(paste(italic("Mpeg1"))))+
  ylim ( 1,6.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(5.75, 4.75),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))

c2p





c1 <- dp.l[[5]]
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
  ggtitle(substitute(paste(italic(" Sell"))))+
  ylim ( -3, 4.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(3.5,2.5),
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
  ggtitle(substitute(paste(italic("Olfr127"))))+
  ylim ( -1,5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(4.25, 3.25),
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
  ggtitle(substitute(paste(italic("Lrrc25"))))+
  ylim ( -1,5.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(4.75, 3.75),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "***"))

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
  ggtitle(substitute(paste(italic("Six6"))))+
  ylim ( -2, 6)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", 'ASC'),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(5, 4),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))

c6p


subdown <- gridExtra::grid.arrange(c1p, c2p, c5p, c4p, c3p, c6p, nrow = 1)

ggsave("manuscript/brain/manuscript70/results/results_figures/ASC_mPFC_DOWN_TOPRegulated_genes_manuscriptSignifcants.png",subdown, height = 3, width= 18,  dpi=300)
