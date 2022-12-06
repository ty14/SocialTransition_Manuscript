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

cdes <- limma_list$controldes 

domdes <- limma_list$domdes 


#looking further into descending genes:
des <- read_csv('manuscript/brain/manuscript70/results/tables/descenders_tranisiton_MEA_genes.csv')
head(des)

des %>% filter(reg == "Down")

ds <- des$symbol

cdes_top <- cdes[cdes$symbol %in% des$symbol,]

cd_down <- cdes_top %>% arrange(-logFC)
head(cd_down, 10)

dd_top <- domdes[domdes$symbol %in% des$symbol,]
dd_down <- dd_top %>% arrange(logFC)
head(dd_down, 10)

sd <- dd_top

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


gettop10GO(des, my_showCategory) -> top10go1

# write.csv(top10go1, 'manuscript/brain/manuscript70/results/tables/Descenders_GOterms_BP.csv' ,row.names = F)


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


up <- top10go1 %>% arrange(qvalue) %>% filter(direction == "Up") 
down <- top10go1%>% arrange(qvalue) %>% filter(direction == "Down")
simMatrix <- calculateSimMatrix(down$ID,
                                orgdb="org.Mm.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10(down$qvalue), down$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=.99,
                                orgdb="org.Mm.eg.db")

scatterPlot(simMatrix, reducedTerms) 



### Make DGE boxplots:

ex <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_MeA_ControlDD.RDS")

dm <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

dom_m_id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

a1<- dm %>% filter(symbol == "Mybpc1")
a2 <- dm %>% filter(symbol == "Cldn3")
a3 <- dm %>% filter(symbol == "Mas1")
# a4 <- dm %>% filter(symbol == "Hcrt")
a5<- dm %>% filter(symbol == "Crym")
a6 <- dm %>% filter(symbol == "Emx1")
a7 <- dm %>% filter(symbol == "Vip")
a8 <- dm %>% filter(symbol == "Spink8")
# a9<- dm %>% filter(symbol == "Hgf")
# a10 <- dm %>% filter(symbol == "Egr2")
a11 <- dm %>% filter(symbol == "Tex15")
a12 <- dm %>% filter(symbol == "Bhlhe22")
a13 <- dm %>% filter(symbol == "Ptgs2")
a14 <- dm %>% filter(symbol == "Tal1")
# a15 <- dm %>% filter(symbol == "Fezf2")

a1<- dm %>% filter(symbol == "Unc13a")
a2 <- dm %>% filter(symbol == "Ddn")
a3 <- dm %>% filter(symbol == "Tspan17")

dp <- a1 %>% rbind(a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15)
dp1 <- a1 %>%  rbind(a2, a3)

colnames(d)

dp2 <- dp1  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(dom_m_id)

dp2$group <- factor(dp2$group, levels = c("CDOM", "DOM", "DES")) 


source('functions/geom_boxjitter.R')

p1 <- ggplot(dp2, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  facet_wrap(~symbol,scales="free_y", ncol = 6)+
  ylab("MEA Normalized Expression") +
  xlab("")+
  theme_bw()+
  theme(legend.position = "none")
p1

ggsave("manuscript/brain/manuscript70/results/results_figures/DES_TOPRegulated_genes_manuscript.png",p1, height = 2.25, width= 12,  dpi=300)


dp.l <- split(dp2 , dp2$symbol)
#Cldn3 c **, d **
#Crym c**, d ***
#Mas1 c **, d **
#Mybpc1 c **, d *
#spink8 c *, **
#vip c *, d * 


c1 <- dp.l[[3]]
 c1p <-  ggplot(c1, aes(group,value, color = group, fill = group))+
    geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                   alpha = 0.4, jitter.size = 2,, size = .8,
                   jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                   position = position_dodge(0.85)) +
    scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
    scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
   ylim(- 8, max(c1$value) +2)+
    ylab("") +
    xlab("")+
   ggtitle(substitute(paste(italic("Unc13a"))))+
   ylim ( 5,9.5)+
    theme_bw(base_size = 13)+
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    geom_signif(comparisons = list(c("CDOM", "DES"),
                                   c("DOM", "DES")),
             map_signif_level = TRUE,
             y_position = c(8.5, 7.5),
             col = 2,
             size = 1,
             textsize = 5,
             annotations = c("**", "**"))

 c1p
 c5p
 c4p
 c3p
 c2p
 c1p
 

 desup <- gridExtra::grid.arrange(c1p, c2p, c3p, c4p, c5p, c6p, nrow = 1)
 
 ggsave("manuscript/brain/manuscript70/results/results_figures/DES_TOPRegulated_genes_manuscriptSignifcants.png",desup, height = 3, width= 18,  dpi=300)
 
 
 ##Down regulated genes
a1 <- dm %>% filter(symbol == "Grb7")
 a2 <- dm %>% filter(symbol == "Slitrk6")
 a3 <- dm %>% filter(symbol == "Erbb3")
 a4 <- dm %>% filter(symbol == "Fa2h")
 a5<- dm %>% filter(symbol == "Dmrtb1")
 a6 <- dm %>% filter(symbol == "Serpinb1a")
 a7 <- dm %>% filter(symbol == "Mog")
 a8 <- dm %>% filter(symbol == "Lpar1")
a9<- dm %>% filter(symbol == "Mobp")
a10 <- dm %>% filter(symbol == "Gjc2")
 a11 <- dm %>% filter(symbol == "Fbln5")
 a12 <- dm %>% filter(symbol == "Mab21l1")
 a13 <- dm %>% filter(symbol == "Col14a1")
 a14 <- dm %>% filter(symbol == "1700066B19Rik")
 a15 <- dm %>% filter(symbol == "Tmem254a")
 
 dp <- a1 %>% rbind(a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15)
 dp1 <- a1 %>%  rbind(a2, a3, a5, a4,a11)
 
 
 dp2 <- dp1  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(dom_m_id)
 
 dp2$group <- factor(dp2$group, levels = c("CDOM", "DOM", "DES")) 
 
 p1 <- ggplot(dp2, aes(group,value, color = group, fill = group))+
   geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                  alpha = 0.4, jitter.size = 2,, size = .8,
                  jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                  position = position_dodge(0.85)) +
   scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
   scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
   facet_wrap(~symbol,scales="free_y", ncol = 6)+
   ylab("MEA Normalized Expression") +
   xlab("")+
   theme_bw()+
   theme(legend.position = "none")
 p1
 
 ggsave("manuscript/brain/manuscript70/results/results_figures/DES_TOP_Down_Regulated_genes_manuscript.png",p1, height = 2.5, width= 12,  dpi=300)

 #Dmrtb1, c*, d **
 #Erbb3, c*, d **
 #Fa2h, c**, d ***
 #fbln5, c**, d *
 #Grb7, c*, d **
 #Slitrk6, c***, d ***
 
 dp.l <- split(dp2 , dp2$symbol)
 
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
   ggtitle(substitute(paste(italic("Slitrk6"))))+
   ylim ( -1,6)+
   theme_bw(base_size = 13)+
   theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
   geom_signif(comparisons = list(c("CDOM", "DES"),
                                  c("DOM", "DES")),
               map_signif_level = TRUE,
               y_position = c(5.5, 4.5),
               col = 2,
               size = 1,
               textsize = 5,
               annotations = c("***", "***"))
 
 c6p
 c5p
 c4p
 c3p
 c2p
 c1p
 
 
 desdown <- gridExtra::grid.arrange(c1p, c2p, c3p, c4p, c5p, c6p, nrow = 1)
 ggsave("manuscript/brain/manuscript70/results/results_figures/DES_TOP_down_Regulated_genes_manuscriptSignifcants.png",desdown, height = 3, width= 18,  dpi=300)
 
 
 
 
 ########## cluster go terms stuff 
 
 df <- read_csv('manuscript/brain/manuscript70/results/tables/Descenders_GOterms_BP.csv')

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
 df_upx.filename <-file.path("manuscript/brain/manuscript70/results/tables/MEA_DOM_up30_clusterprofiler_cluster_enr_results.txt")
 write.table(df_upx,df_upx.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)

 
  
 colnames(df_downx ) <- c("Name","Description", "pvalue","qvalue", "phenotype", "genes")
 head(df_downx)
 # write file for cytoscape
 df_downx.filename <-file.path("manuscript/brain/manuscript70/results/tables/MEA_DOM_down30_clusterprofiler_cluster_enr_results.txt")
 write.table(df_downx,df_downx.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)
 
 
 
 ########### heat plot
 
 des <- read_csv('manuscript/brain/manuscript70/results/tables/descenders_tranisiton_mPFC_genes.csv')
 head(des)
 
 sd <- domdes[domdes$symbol %in% des$symbol,]
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






