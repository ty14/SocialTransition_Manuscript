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

limma_list<- readRDS("brain/results/RDS/limma_MEA_CSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

y1a <- limma_list$csub

y2a <- limma_list$casc

y3a <- limma_list$subasc


limma_list<- readRDS("brain/results/RDS/limma_mPFC_CSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 


y1b <- limma_list$csub

y2b <- limma_list$casc
 
y3b <- limma_list$subasc

my_ont = "BP"
my_showCategory = 100

gettop10GO <- function(limma_df,my_showCategory = 10){
  go_df_up <- limma_df %>% 
    filter(logFC>0) %>% 
    left_join(grcm38 %>% dplyr::select(symbol, entrez)) %>% 
    mutate(entrez = as.character(entrez)) %>% 
    dplyr::select(entrez, logFC) %>% 
    filter(!is.na(entrez)) %>% 
    arrange(desc(logFC))
  
  
  go_df_down <- limma_df %>% 
    filter(logFC<0) %>% 
    left_join(grcm38 %>% dplyr::select(symbol, entrez)) %>% 
    mutate(entrez = as.character(entrez)) %>% 
    dplyr::select(entrez, logFC) %>% 
    filter(!is.na(entrez)) %>% 
    arrange(desc(logFC))
  
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
                       pvalueCutoff  = 0.2, #
                       qvalueCutoff  = 0.20) 
  
  fortify(
    ggo_up,
    showCategory = my_showCategory,
    by = "GeneRatio",
    split = NULL,
    includeAll = TRUE
  ) %>% 
    arrange(desc(GeneRatio)) %>% 
    mutate(direction = "Up") -> temp1
  
  
  fortify(
    ggo_down,
    showCategory = my_showCategory,
    by = "GeneRatio",
    split = NULL,
    includeAll = TRUE
  ) %>% 
    arrange(desc(GeneRatio)) %>% 
    mutate(direction = "Down") -> temp2
  return(rbind(temp1,temp2))
  
}


gettop10GO(y1a, my_showCategory) %>% 
  mutate(comparison = "CSUB-SUB") -> top10go1

gettop10GO(y2a, my_showCategory ) %>% 
  mutate(comparison = "CSUB-ASC") -> top10go2

gettop10GO(y3a, my_showCategory ) %>% 
  mutate(comparison = "SUB-ASC") -> top10go3


rbind(top10go1,top10go2,top10go3) -> top10_GOterms

write.csv(top10_GOterms,"brain/results/tables/topBP_GOterms_MeA_CSUB.csv", row.names = F)


gettop10GO(y1b, my_showCategory) %>% 
  mutate(comparison = "CSUB-SUB") -> top10go1

gettop10GO(y2b, my_showCategory ) %>% 
  mutate(comparison = "CSUB-ASC") -> top10go2

gettop10GO(y3b, my_showCategory ) %>% 
  mutate(comparison = "SUB-ASC") -> top10go3

rbind(top10go1,top10go2,top10go3) -> top10_GOterms

 write.csv(top10_GOterms,"brain/results/tables/topBP_GOterms_mPFC_CSUB.csv", row.names = F)



up <- top10go1 %>% arrange(qvalue)  %>%  head(20) %>% filter(direction == "Up")
# view

top10_GOterms$Description[1:5]
top10_GOterms$Description[6:10]
top10_GOterms$Description[11:15]
top10_GOterms$Description[16:20]
top10_GOterms$Description[21:25]


