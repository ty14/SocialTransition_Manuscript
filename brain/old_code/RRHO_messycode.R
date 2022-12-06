#libraries 
library(tidyverse)
library(RRHO2)



my_logFC_threshold = 0.2

#MEA data

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_MEA_ControlDD.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol,ensgene, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

cdom <- limma_list$controldom

cdes <- limma_list$controldes 

domdes <- limma_list$domdes 


#looking further into descending genes:
des <- read_csv('manuscript/brain/manuscript70/results/tables/descenders_tranisiton_MEA_genes.csv')
head(des)


sd <- domdes[domdes$symbol %in% des$symbol,]


#MEA data ## sub 

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_MEA_ControlSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

cs <- limma_list$controlsub

ca <- limma_list$controlasc 

sa <- limma_list$subasc 


#looking further into descending genes:
asc <- read_csv('manuscript/brain/manuscript70/results/tables/ascenders_tranisiton_MEA_genes.csv')
head(asc)


asc <- sa[sa$symbol %in% asc$symbol,]

####################

sdx <- cdom[cdom$symbol %in% cs$symbol,]

sdx <- sdx[!duplicated(sdx[ , c("symbol")]),]

ascx <- cs[cs$symbol%in% cdom$symbol,]

ascx  <- ascx[!duplicated(ascx[ , c("symbol")]),]

SD <- sdx %>%  mutate(gene = row_number())
AS <- ascx %>% mutate(gene = row_number())




SD$gene <- paste0("Gene", SD$gene)
head(SD)
AS$gene <- paste0("Gene", AS$gene)
head(AS)
#List 1
Gene = SD$gene

list1_pvalue_1_200 <- SD%>% filter(logFC > 0 & P.Value < 0.01) 
list1_pvalue_1_200 <-list1_pvalue_1_200$P.Value## up-regulated genes
list1_pvalue_201_400 <- SD%>% filter(logFC < 0 & P.Value < 0.01)
list1_pvalue_201_400 <-list1_pvalue_201_400$P.Value ## down-regulated genes
an <- subset(SD, !(P.Value %in% list1_pvalue_1_200))
an <- subset(an, !(P.Value %in% list1_pvalue_201_400))
list1_pvalue_401_2000 <- an$P.Value  ## non-changed genes
list1_DDE <- c(-log10(list1_pvalue_1_200), 
               -log10(list1_pvalue_201_400) * (-1), 
               -log10(list1_pvalue_401_2000) * sample(c(1,-1), length(list1_pvalue_401_2000), replace = TRUE))

lapply(list1_DDE, head)
gene_list1 <- data.frame(Genes=Gene,DDE = list1_DDE, stringsAsFactors = FALSE)

#List2
list2_pvalue_1_200 <- AS %>% filter(logFC > 0 & P.Value < 0.01) 
list2_pvalue_1_200 <-list2_pvalue_1_200$P.Value## up-regulated genes
list2_pvalue_201_400 <- AS%>% filter(logFC < 0 & P.Value < 0.01)
list2_pvalue_201_400 <-list2_pvalue_201_400$P.Value ## down-regulated genes
pn <- subset(AS, !(P.Value %in% list2_pvalue_1_200))
pn <- subset(pn, !(P.Value %in% list2_pvalue_201_400))
list2_pvalue_401_2000 <- pn$P.Value  ## non-changed genes
list2_DDE <- c(-log10(list2_pvalue_1_200), 
               -log10(list2_pvalue_201_400) * (-1), 
               -log10(list2_pvalue_401_2000) * sample(c(1,-1), length(list2_pvalue_401_2000), replace = TRUE))


gene_list2 <- data.frame(Genes=Gene,DDE = list2_DDE, stringsAsFactors = FALSE)

library(RRHO2)

RRHO_obj <-  RRHO2_initialize(gene_list1, gene_list2, labels = c("CDOM vs. DOM", "CSUB vs. SUB"), log10.ind=TRUE)


RRHO2_heatmap(RRHO_obj)

