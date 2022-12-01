library(Mus.musculus)
organism = 'org.Mm.eg.db'
library(edgeR)
library(clusterProfiler)
library(gprofiler2)
library(tidyverse)


my_logFC_threshold = 0.2

#MEA data

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_MEA_ControlDD.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol,ensgene, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

domdes <- limma_list$domdes


des <- read_csv('manuscript/brain/manuscript70/results/tables/descenders_tranisiton_MEA_genes.csv')
head(des)

sd <- domdes[domdes$symbol %in% des$symbol,]
head(sd)

sd <- sd[!duplicated(sd[ , c("symbol")]),]

# write.csv(sd, 'manuscript/brain/manuscript70/results/tables/sd_results_forEnrichMAP.csv' ,row.names = F)
# write.table(sd,'manuscript/brain/manuscript70/results/tables/sd_results_forEnrichMAP.txt',col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)

sd_upreg <- sd %>% filter(logFC <= -0.2)

sd_downreg <- sd %>% filter(logFC >= 0.2)



# # we want the log2 fold change 
original_gene_list_up <- sd_upreg$logFC

original_gene_list_down <- sd_downreg$logFC

original_gene_list <- sd$logFC
# name the vector
names(original_gene_list_up) <- sd_upreg$entrez

names(original_gene_list_down) <- sd_downreg$entrez


names(original_gene_list) <- sd$ensgene
# omit any NA values
up_list<-na.omit(original_gene_list_up)

down_list<-na.omit(original_gene_list_down)

gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
up_list = sort(up_list, decreasing = TRUE)

down_list = sort(down_list, decreasing = TRUE)

gene_list = sort(gene_list, decreasing = TRUE)

up<- enrichGO(gene = sd_upreg$entrez %>% unique(),
                  OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  readable = T,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.2,
                  qvalueCutoff  = 0.20)


down<- enrichGO(gene = sd_downreg$entrez %>% unique(),
              OrgDb = org.Mm.eg.db::org.Mm.eg.db,
              keyType = "ENTREZID",
              ont = "BP",
              readable = T,
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.2,
              qvalueCutoff  = 0.20)



library(gprofiler2)



gostres <- gost(query = des_up, 
                organism = "mmusculus", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = "GO:BP", as_short_link = FALSE)


sd_up <- gostres[[1]] %>% as.data.frame()

sd_upfilename <-file.path("manuscript/brain/manuscript70/results/tables/MEA_DOM_up_gPROFILER_cluster_enr_results.txt")
write.table(sd_up,sd_upfilename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)

sd_down <- gostres[[1]]
sd_downfilename <-file.path("manuscript/brain/manuscript70/results/tables/MEA_DOM_down_gPROFILER_cluster_enr_results.txt")
write.table(sd_down,sd_downfilename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)



## get results 
# egobp.results.df <- up@result %>% as.data.frame()
egobp.results.df <- down@result %>% as.data.frame()
# create a new column for term size from BgRatio
egobp.results.df$term.size <- gsub('/(\\d+)*', "", egobp.results.df$BgRatio)

# filter for term size to keep only term.size => 3, gene count >= 5 and subset
egobp.results.df <- egobp.results.df[which(egobp.results.df[,'term.size'] >= 3 & egobp.results.df[,'Count'] >= 5),]
egobp.results.df <- egobp.results.df[c("ID", "Description", "pvalue", "qvalue", "geneID")]

# format gene list column
egobp.results.df$geneID <- gsub("/", ",", egobp.results.df$geneID)

# add column for phenotype
egobp.results.df <- cbind(egobp.results.df, phenotype= -1)
egobp.results.df <- egobp.results.df[, c(1, 2, 3, 4, 6, 5)]

# change column headers
colnames(egobp.results.df) <- c("Name","Description", "pvalue","qvalue", "phenotype", "genes")
# write_csv <-file.path(getwd(),paste("clusterprofiler_cluster_enr_results.txt",sep="_"))
write.csv(egobp.results.df, 'manuscript/brain/manuscript70/results/tables/MEA_DOM_down_clusterprofiler_cluster_enr_results.csv' ,row.names = F)


egobp.results.filename <-file.path("manuscript/brain/manuscript70/results/tables/MEA_DOM_down_clusterprofiler_cluster_enr_results.txt")
write.table(egobp.results.df,egobp.results.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)



gse <- gseGO(geneList=gene_list, 
             ont ="BP",
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")


gse2 <- pairwise_termsim(gse)
set.seed(123)
p1 <- emapplot(gse2,color = "enrichmentScore", layout="kk", group_category = T, group_legend = F, 
               cex_label_group = 1.5, node_label = "group", nCluster = 10) 

p1

set.seed(123)
emapplot(
  gse2,
  showCategory = 30,
  layout = NULL,
  coords = NULL,
  color = "enrichmentScore",
  min_edge = 0.6,
  cex_label_category = 1,
  cex_category = 1,
  cex_line = 1,
  shadowtext = TRUE,
  label_style = "ggforce",
  repel = FALSE,
  node_label = "group",
  with_edge = TRUE,
  group_category =T,
  group_legend = FALSE,
  cex_label_group = 1,
  nWords = 4,
  label_format = 16,
  clusterFunction = stats::kmeans,
  nCluster = 5
)




cnetplot(gse2, node_label="category", 
         cex_label_category = 1.2)


heatplot(gse2, foldChange=gene_list, showCategory=10)
