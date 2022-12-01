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

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_MEA_ControlDD.RDS") %>% 
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
colnames(up_des)<-"symbol"

up_dom<- up_dom %>% mutate(reg = "Up")

cddd <- cdom %>% filter(reg == "up")

dddd <- domdes %>% filter(reg == "down")

down_dom <- cddd$symbol[(cddd$symbol %in% dddd$symbol)] %>% as.data.frame(.)
colnames(down_des)<-"symbol"

down_dom <- down_dom %>% mutate(reg = "Down")

#up dom = 451, down dom = 375
dom <- up_dom %>% full_join(down_dom)

dom$symbol <- dom$.

dom <- dom %>% select(-.)

write.csv(dom, 'manuscript/brain/manuscript70/results/tables/dominant_distruption_MEA_genes.csv' ,row.names = F)


#looking further into descending genes:
dom <- read_csv('manuscript/brain/manuscript70/results/tables/dominant_distruption_MEA_genes.csv')
head(dom)



cdom_top <- cdom[cdom$symbol %in% dom$symbol,]

cd_down <- cdom_top %>% arrange(-logFC)
x <- head(cd_down, 20)

dd_top <- domdes[domdes$symbol %in% dom$symbol,]
dd_down <- dd_top %>% arrange(logFC)
y <-head(dd_down, 20)
x[x$symbol %in% y$symbol,]


# UP: Pvalb, lhx8,Slc18a3, Esrrb, Slc10a4, Chrm2,

#DOWN: Col5a1, Esyt3, Trh, Smpdl3b, Ankef1, Krt77


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

write.csv(top10go1, 'manuscript/brain/manuscript70/results/tables/Dominant_GOterms_BP.csv' ,row.names = F)


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

ex <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_MeA_ControlDD.RDS")

sm <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

sub_m_id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

# UP: Pvalb, lhx8,Slc18a3, Esrrb, Slc10a4, Chrm2,

a1<- sm %>% filter(symbol == "Pvalb")
a2 <- sm %>% filter(symbol == "Lhx8")
a3 <- sm %>% filter(symbol == "Slc18a3")
a4 <- sm %>% filter(symbol == "Esrrb")
a5<- sm %>% filter(symbol == "Slc10a4")
a6 <- sm %>% filter(symbol == "Chrm2")


dp <- a1 %>% rbind(a2,a3,a4,a5,a6)


colnames(dp)

dp2 <- dp  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(sub_m_id)

dp2$group <- factor(dp2$group, levels = c("CDOM", "DOM", "DES")) 


source('functions/geom_boxjitter.R')


dp.l <- split(dp2 , dp2$symbol)
#Pvalb c ***, d ***
#Lhx8 c*, d **
#Slc18a3 c *, d **
#Esrrb c *, d *
#Slc10a4 c *, ***
#Chrm2 c **, d *** 
# MEA Normalized Expression

c1 <- dp.l[[1]]
c1p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("MEA Normalized Expression") +
  xlab("")+
  ggtitle(substitute(paste(italic("Chrm2"))))+
  ylim ( 2,8)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(7.5, 6.5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("**", "***"))

c1p



c1 <- dp.l[[2]]
c2p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Esrrb"))))+
  ylim ( -3,6.5)+
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

c2p



c1 <- dp.l[[3]]
c3p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Lhx8"))))+
  ylim ( -3,8.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(7.5, 6.5),
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
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Slc10a4"))))+
  ylim ( 0,8)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(7.25, 6.5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "***"))

c4p


c1 <- dp.l[[6]]
c5p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Slc18a3"))))+
  ylim ( -3,6)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(5.25, 4.25),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "**"))

c5p


c1 <- dp.l[[4]]
c6p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Pvalb"))))+
  ylim ( 1,9.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(8.75,8),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("***", "***"))

c6p


domup <- gridExtra::grid.arrange(c1p, c2p, c5p, c4p, c3p, c6p, nrow = 1)

ggsave("manuscript/brain/manuscript70/results/results_figures/DOM_TOPRegulated_genes_manuscriptSignifcants.png",domup, height = 3, width= 18,  dpi=300)

# DOWN: Col5a1, Esyt3, Trh, Smpdl3b, Ankef1, Krt77


a1<- sm %>% filter(symbol == "Col5a1")
a2 <- sm %>% filter(symbol == "Esyt3")
a3 <- sm %>% filter(symbol == "Trh")
a4 <- sm %>% filter(symbol == "Smpdl3b")
a5<- sm %>% filter(symbol == "Ankef1")
a6 <- sm %>% filter(symbol == "Krt77")


dp <- a1 %>% rbind(a2,a3,a4,a5,a6)


colnames(dp)

dp2 <- dp  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(sub_m_id)

dp2$group <- factor(dp2$group, levels = c("CDOM", "DOM", "DES")) 


source('functions/geom_boxjitter.R')


dp.l <- split(dp2 , dp2$symbol)
# Col5a1, c ***, des**
# Esyt3,  **, **
# Trh, *, **
# Smpdl3b, *,**
# Ankef1, *, *
# Krt77 ***, ***


c1 <- dp.l[[1]]
c1p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab(" MEA Normalized Expression") +
  xlab("")+
  ggtitle(substitute(paste(italic("Ankef1"))))+
  ylim ( -2 ,6)+
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
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Col5a1"))))+
  ylim ( -2.5,5.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(4.5, 3.5),
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
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Esyt3"))))+
  ylim ( -1,7)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(6, 5),
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
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Krt77"))))+
  ylim ( 2,7.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(6.75,5.75),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("***", "***"))


c4p



c1 <- dp.l[[5]]
c5p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Smpdl3b"))))+
  ylim ( -1,7.5)+
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


c5p



c1 <- dp.l[[6]]
c6p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("") +
  xlab("")+
  ggtitle(substitute(paste(italic("Trh"))))+
  ylim ( 1,8.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CDOM", "DOM"),
                                 c("DOM", "DES")),
              map_signif_level = TRUE,
              y_position = c(7.5, 6.5),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "**"))


c6p



domdown <- gridExtra::grid.arrange(c1p, c2p, c3p, c4p, c5p, c6p, nrow = 1)




ggsave("manuscript/brain/manuscript70/results/results_figures/DOM_DOWNRegulated_genes_manuscriptSignifcants.png",domdown, height = 3, width= 18,  dpi=300)



########## cluster go terms stuff 

df <- read_csv('manuscript/brain/manuscript70/results/tables/Dominant_GOterms_BP.csv')

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
df_upx.filename <-file.path("manuscript/brain/manuscript70/results/tables/MEA_DDOM_up30_clusterprofiler_cluster_enr_results.txt")
write.table(df_upx,df_upx.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)



colnames(df_downx ) <- c("Name","Description", "pvalue","qvalue", "phenotype", "genes")
head(df_downx)
# write file for cytoscape
df_downx.filename <-file.path("manuscript/brain/manuscript70/results/tables/MEA_DDOM_down30_clusterprofiler_cluster_enr_results.txt")
write.table(df_downx,df_downx.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)



########### heat plot

dom <- read_csv('manuscript/brain/manuscript70/results/tables/dominant_distruption_MEA_genes.csv')
head(dom)

sd <- domdes[cdom$symbol %in% dom$symbol,]
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


