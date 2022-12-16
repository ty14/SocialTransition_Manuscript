# libraries 
library(limma)
library(edgeR)
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
library(tidyverse)
grcm38 # mouse genes


my_logFC_threshold = 0.2

y1a <- read_csv("brain/results/tables/descender_MEA_genes.csv")

y1b <- read_csv("brain/results/tables/descender_mPFC_genes.csv")


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



gettop10GO(y1a, my_showCategory) %>% 
  mutate(comparison = "DES - MEA") -> top10go1

gettop10GO(y1b, my_showCategory ) %>% 
  mutate(comparison = "DES - mPFC") -> top10go2



rbind(top10go1,top10go2) -> top10_GOterms

write.csv(top10_GOterms,"brain/results/tables/topBP_GOterms_DES_both.csv", row.names = F)

# view

top10_GOterms$Description[1:5]
top10_GOterms$Description[6:10]
top10_GOterms$Description[11:15]
top10_GOterms$Description[16:20]
top10_GOterms$Description[21:25]

######### Ascenders

y2a <- read_csv("brain/results/tables/ascender_MEA_genes.csv")

y2b <- read_csv("brain/results/tables/ascender_mPFC_genes.csv")


gettop10GO(y1a, my_showCategory) %>% 
  mutate(comparison = "ASC - MEA") -> top10go1

gettop10GO(y2b, my_showCategory ) %>% 
  mutate(comparison = "ASC - mPFC") -> top10go2



rbind(top10go1,top10go2) -> top10_GOterms

write.csv(top10_GOterms,"brain/results/tables/topBP_GOterms_ASC_both.csv", row.names = F)

# view

top10_GOterms$Description[1:5]
top10_GOterms$Description[6:10]
top10_GOterms$Description[11:15]
top10_GOterms$Description[16:20]
top10_GOterms$Description[21:25]


###cdom

y2a <- read_csv("brain/results/tables/distruption_MEA_genes.csv")

y2b <- read_csv("brain/results/tables/distruption_mPFC_genes.csv")


gettop10GO(y2a, my_showCategory) %>% 
  mutate(comparison = "DOM DIS - MEA") -> top10go1

gettop10GO(y2b, my_showCategory ) %>% 
  mutate(comparison = " DOM DIS - mPFC") -> top10go2



rbind(top10go1,top10go2) -> top10_GOterms

write.csv(top10_GOterms,"brain/results/tables/topBP_GOterms_DOM_DIS_both.csv", row.names = F)

# view

top10_GOterms$Description[1:5]
top10_GOterms$Description[6:10]
top10_GOterms$Description[11:15]
top10_GOterms$Description[16:20]
top10_GOterms$Description[21:25]


###cdom

y2a <- read_csv("brain/results/tables/distruptionSUB_MEA_genes.csv")

y2b <- read_csv("brain/results/tables/distruptionSUB_mPFC_genes.csv")


gettop10GO(y2a, my_showCategory) %>% 
  mutate(comparison = "SUB DIS - MEA") -> top10go1

gettop10GO(y2b, my_showCategory ) %>% 
  mutate(comparison = " SUB DIS - mPFC") -> top10go2



rbind(top10go1,top10go2) -> top10_GOterms

write.csv(top10_GOterms,"brain/results/tables/topBP_GOterms_SUB_DIS_both.csv", row.names = F)

# view

top10_GOterms$Description[1:5]
top10_GOterms$Description[6:10]
top10_GOterms$Description[11:15]
top10_GOterms$Description[16:20]
top10_GOterms$Description[21:25]


