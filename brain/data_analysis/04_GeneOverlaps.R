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


my_logFC_threshold = 0.2

limma_list<- readRDS("brain/results/RDS/limma_MEA_CDOM.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez)))


cdom <- limma_list$cdom
cdom$reg <- ifelse(cdom$logFC >= 0.2, "up", "down")

cdes <- limma_list$cdes 
cdes$reg <- ifelse(cdes$logFC >= 0.2, "up", "down")

domdes <- limma_list$domdes 
domdes$reg <- ifelse(domdes$logFC >= 0.2, "up", "down")

#looking for des upregulated genes: 
cdd <- cdes %>% filter(reg == "down")

ddd <- domdes %>% filter(reg == "down")

up_des <- cdd$symbol[(cdd$symbol %in% ddd$symbol)] %>% as.data.frame(.)
colnames(up_des)<-"symbol"

up_des <- up_des %>% mutate(reg = "Up")

#looking for des downregulated genes:
cddd <- cdes %>% filter(reg == "up")

dddd <- domdes %>% filter(reg == "up")

down_des <- cddd$symbol[(cddd$symbol %in% dddd$symbol)] %>% as.data.frame(.)
colnames(down_des)<-"symbol"

down_des <- down_des %>% mutate(reg = "Down")

#up descenders is 206, down decenders is 160 #MEA
#up descenders is 55, down decenders is 56 #mPFC
des <- up_des %>% full_join(down_des)

#save
# write.csv(des, 'brain/results/tables/descender_MEA_genes.csv' ,row.names = F)
# write.csv(des, 'brain/results/tables/descender_mPFC_genes.csv' ,row.names = F)


############ 

## looking for genes different from CDOM  

# downreg in distrupted animals : 
cdd <- cdes %>% filter(reg == "up")

cddom <- cdom %>% filter(reg == "up")

down_cdom <- cdd$symbol[(cdd$symbol %in% cddom$symbol)] %>% as.data.frame(.)
colnames(down_cdom)<-"symbol"

down_dis <- down_cdom %>% mutate(reg = "Down")

# upreg in distrupted animals :
cddd <- cdes %>% filter(reg == "down")

cdddom <- cdom %>% filter(reg == "down")

up_cdom <- cddd$symbol[(cddd$symbol %in% cdddom$symbol)] %>% as.data.frame(.)
colnames(up_cdom)<-"symbol"

up_dis <- up_cdom %>% mutate(reg = "Up")

#up cdom is 115, down cdom is 98
cdom_gene <- up_cdom %>% full_join(down_cdom)


#up dis is 93, down decenders is 103
#up dis is 96, down decenders is 94
dis <- up_dis %>% full_join(down_dis)

#save
# write.csv(dis, 'brain/results/tables/distruption_MEA_genes.csv' ,row.names = F)
# write.csv(dis, 'brain/results/tables/distruption_mPFC_genes.csv' ,row.names = F)

##Tmem42 is only gene overlaps between dis & des??? # MEA ; none in mPFC 
dis$symbol[(dis$symbol %in% des$symbol)]

############################################
#Ascenders 

my_logFC_threshold = 0.2

limma_list<- readRDS("brain/results/RDS/limma_mPFC_CSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez)))


csub <- limma_list$csub
csub$reg <- ifelse(csub$logFC >= 0.2, "up", "down")

casc <- limma_list$casc 
casc$reg <- ifelse(casc$logFC >= 0.2, "up", "down")

subasc <- limma_list$subasc 
subasc$reg <- ifelse(subasc$logFC >= 0.2, "up", "down")

#looking for asc upregulated genes: 
cdd <- casc %>% filter(reg == "down")

ddd <- subasc %>% filter(reg == "down")

up_asc <- cdd$symbol[(cdd$symbol %in% ddd$symbol)] %>% as.data.frame(.)
colnames(up_asc)<-"symbol"

up_asc <- up_asc %>% mutate(reg = "Up")

#looking for asc downregulated genes:
cddd <- casc %>% filter(reg == "up")

dddd <- subasc %>% filter(reg == "up")

down_asc <- cddd$symbol[(cddd$symbol %in% dddd$symbol)] %>% as.data.frame(.)
colnames(down_asc)<-"symbol"

down_asc <- down_asc %>% mutate(reg = "Down")

#up ascenders is 88, down ascenders is 67 #MEA
#up ascenders is 56, down ascenders is 64 #mPFC
asc <- up_asc %>% full_join(down_asc)

#save
# write.csv(asc, 'brain/results/tables/ascender_MEA_genes.csv' ,row.names = F)
# write.csv(asc, 'brain/results/tables/ascender_mPFC_genes.csv' ,row.names = F)


############ 

## looking for genes different from CSUB 

# downreg in distrupted animals : 
cdd <- casc %>% filter(reg == "up")

cdsub <- csub %>% filter(reg == "up")

down_csub <- cdd$symbol[(cdd$symbol %in% cdsub$symbol)] %>% as.data.frame(.)
colnames(down_csub)<-"symbol"

down_dis <- down_csub %>% mutate(reg = "Down")

# upreg in distrupted animals :
cddd <- casc %>% filter(reg == "down")

cddsub <- csub %>% filter(reg == "down")

up_csub <- cddd$symbol[(cddd$symbol %in% cddsub$symbol)] %>% as.data.frame(.)
colnames(up_csub)<-"symbol"

up_dis <- up_csub %>% mutate(reg = "Up")


#up dis is 39, down decenders is 31 #MEA
#up dis is 101, down decenders is 110 #mPFC
dis <- up_dis %>% full_join(down_dis)

#save
# write.csv(dis, 'brain/results/tables/distruptionSUB_MEA_genes.csv' ,row.names = F)
# write.csv(dis, 'brain/results/tables/distruptionSUB_mPFC_genes.csv' ,row.names = F)

##Lrba is only gene overlaps between dis & asc??? # MEA ; none in mPFC 
asc[(asc$symbol %in% dis$symbol)]

