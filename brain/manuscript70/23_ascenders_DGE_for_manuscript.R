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

ca <- limma_list$controlasc 

sa <- limma_list$subasc 


#looking further into descending genes:
asc <- read_csv('manuscript/brain/manuscript70/results/tables/ascenders_tranisiton_MEA_genes.csv')
head(asc)


casc_top <- ca[ca$symbol %in% asc$symbol,]

cd_up <- casc_top %>% arrange(logFC)
 x <- head(cd_up, 15)

sa_top <- sa[sa$symbol %in% asc$symbol,]
sa_up <- sa_top %>% arrange(logFC)
y <- head(sa_up, 15)


x$symbol[x$symbol %in% y$symbol]


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


gettop10GO(asc, my_showCategory) -> top10go1

write.csv(top10go1, 'manuscript/brain/manuscript70/results/tables/ASC_GOterms_BP.csv' ,row.names = F)

#looking at Transition genes 
# gettop10GO(trans_genes, my_showCategory) -> top10go1 

top10go1$Description[1:5]
top10go1$Description[6:10]
top10go1$Description[11:15]
top10go1$Description[16:20]
top10go1$Description[21:25]

up <- top10go1 %>% arrange(qvalue) %>% filter(direction == "Up")
head(up$Description)

# another way to look at this: 

go_df_up <-  des %>% filter(reg == "Down") %>% 
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

edox <- setReadable(ggo_up, 'org.Mm.eg.db', 'ENTREZID')

p1 <- heatplot(edox, showCategory=5)
p1

# ggsave("manuscript/brain/manuscript70/results/results_figures/DES_GOterms_DownRegulated.png",p1, height = 4, width= 12,  dpi=300)
 
 
library(rrvgo)
head(top10go1)
up <- top10go1 %>% arrange(qvalue) %>% filter(direction == "Up") 
down <- top10go1%>% arrange(qvalue) %>% filter(direction == "Down")
simMatrix <- calculateSimMatrix(up$ID,
                                orgdb="org.Mm.eg.db",
                                ont="BP",
                                method="Rel")

# scores <- setNames(-log10(up$qvalue), up$ID)
# reducedTerms <- reduceSimMatrix(simMatrix,
#                                 scores,
#                                 threshold=.99,
#                                 orgdb="org.Mm.eg.db")
# 
# scatterPlot(simMatrix, reducedTerms) 

### Make DGE boxplots:

ex <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_MeA_ControlSUB.RDS")

sm <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

sub_m_id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

# Tnfrsf25, Ano2, Car12, Adtrp, Vipr1, Adamtsl1, Ikzf1, Acvr1c, Foxo1

a1<- sm %>% filter(symbol == "Unc13a")
a2 <- sm %>% filter(symbol == "Ddn")
a3 <- sm %>% filter(symbol == "Tspan17")
a4 <- sm %>% filter(symbol == "Vipr1")
a5<- sm %>% filter(symbol == "Adamtsl1")
a6 <- sm %>% filter(symbol == "Foxo1")


dp <- a1 %>% rbind(a2,a3)


colnames(dp)

dp2 <- dp  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(sub_m_id)

dp2$group <- factor(dp2$group, levels = c("CSUB", "SUB", "ASC")) 


source('functions/geom_boxjitter.R')


dp.l <- split(dp2 , dp2$symbol)
#Tnfrsf25 c *, d *
#Ano2 c*, d *
#Car12 c *, d **
#Vipr1 c *, d *
#Adamtsl1 c *, *
#foxo1 c *, d ** 
# MEA Normalized Expression

c1 <- dp.l[[2]]
c4p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  ylim(- 8, max(c1$value) +2)+
  ylab("MEA Normalized Expression") +
  xlab("")+
  ggtitle(substitute(paste(italic("Tspan17"))))+
  ylim (4, 8.5)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(8, 7),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "**"))

c6p
c5p
c4p
c3p
c2p
c1p


ascup <- gridExtra::grid.arrange(c1p, c2p, c3p, c4p, c5p, c6p, nrow = 1)

ggsave("manuscript/brain/manuscript70/results/results_figures/ASC_TOPRegulated_genes_manuscriptSignifcants.png",ascup, height = 3, width= 18,  dpi=300)

# Ggt1, Slc6a3, Ngf, Sp8, Tagln, Crhr2, Scx, Iqcd

a1<- sm %>% filter(symbol == "Ggt1")
a2 <- sm %>% filter(symbol == "Crhr2")
a3 <- sm %>% filter(symbol == "Slc6a3")
a4 <- sm %>% filter(symbol == "Sp8")
a5<- sm %>% filter(symbol == "Scx")
a6 <- sm %>% filter(symbol == "Tagln")


dp <- a1 %>% rbind(a2,a3,a4,a5,a6)


colnames(dp)

dp2 <- dp  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(sub_m_id)

dp2$group <- factor(dp2$group, levels = c("CSUB", "SUB", "ASC")) 


source('functions/geom_boxjitter.R')


dp.l <- split(dp2 , dp2$symbol)
# Ggt1, c **, s **
# Crhr2, c * , s*
# Scx, c * , s*
# Slc6a3, c * , s*
# Sp8,  c * , s*
# Tagln,c * , s*


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
  ggtitle(substitute(paste(italic("Tagln"))))+
  ylim ( -1,7)+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("CSUB", "ASC"),
                                 c("SUB", "ASC")),
              map_signif_level = TRUE,
              y_position = c(6.25, 5.25),
              col = 2,
              size = 1,
              textsize = 5,
              annotations = c("*", "*"))

c6p
c5p
c4p
c3p
c2p
c1p


ascdown <- gridExtra::grid.arrange(c1p, c2p, c3p, c4p, c5p, c6p, nrow = 1)




ggsave("manuscript/brain/manuscript70/results/results_figures/ASC_DOWNRegulated_genes_manuscriptSignifcants.png",ascdown, height = 3, width= 18,  dpi=300)



########## cluster go terms stuff 

df <- read_csv('manuscript/brain/manuscript70/results/tables/ASC_GOterms_BP.csv')

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
df_upx.filename <-file.path("manuscript/brain/manuscript70/results/tables/MEA_ASC_up30_clusterprofiler_cluster_enr_results.txt")
write.table(df_upx,df_upx.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)



colnames(df_downx ) <- c("Name","Description", "pvalue","qvalue", "phenotype", "genes")
head(df_downx)
# write file for cytoscape
df_downx.filename <-file.path("manuscript/brain/manuscript70/results/tables/MEA_ASC_down30_clusterprofiler_cluster_enr_results.txt")
write.table(df_downx,df_downx.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)



########### heat plot

asc <- read_csv("manuscript/brain/manuscript70/results/tables/subordinate_distruption_mPFC_genes.csv")
head(asc)

sd <- sa[sa$symbol %in% asc$symbol,]
head(sd)

sd <- sd[!duplicated(sd[ , c("symbol")]),]

# write.csv(sd, 'manuscript/brain/manuscript70/results/tables/sd_results_forEnrichMAP.csv' ,row.names = F)
# write.table(sd,'manuscript/brain/manuscript70/results/tables/sd_results_forEnrichMAP.txt',col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)

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


