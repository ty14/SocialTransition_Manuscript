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

#MEA data

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_MEA_ControlSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

csub <- limma_list$controlsub

csub$reg <- ifelse(csub$logFC >= 0.2, "up", "down")

casc <- limma_list$controlasc
casc$reg <- ifelse(casc$logFC >= 0.2, "up", "down")

sa <- limma_list$subasc
sa$reg <- ifelse(sa$logFC >= 0.2, "up", "down")



#looking for asc upregulated genes: 
cad <- casc %>% filter(reg == "down")

sad <- sa %>% filter(reg == "down")

up_asc <- cad$symbol[(cad$symbol %in% sad$symbol)] %>% as.data.frame(.)
colnames(up_asc)<-"symbol"

up_asc <- up_asc %>% mutate(reg = "Up")

cadd <- casc %>% filter(reg == "up")

sadd <- sa %>% filter(reg == "up")

down_asc <- cadd$symbol[(cadd$symbol %in% sadd$symbol)] %>% as.data.frame(.)
colnames(down_asc)<-"symbol"

down_asc <- down_asc %>% mutate(reg = "Down")

#up asc = 103, down asc = 80
asc <- up_asc %>% full_join(down_asc)


write.csv(asc, 'manuscript/brain/manuscript70/results/tables/ascenders_tranisiton_MEA_genes.csv' ,row.names = F)

## looking for genes different from CSUB: 
cdd <- casc %>% filter(reg == "up")

cds <- csub %>% filter(reg == "up")

up_csub <- cdd$symbol[(cdd$symbol %in% cds$symbol)] %>% as.data.frame(.)
colnames(up_csub)<-"symbol"

up_csub <- up_csub %>% mutate(reg = "Up")

cadd <- casc %>% filter(reg == "down")

cdsub<- csub %>% filter(reg == "down")

down_csub <- cadd$symbol[(cadd$symbol %in% cdsub$symbol)] %>% as.data.frame(.)
colnames(down_csub)<-"symbol"

down_csub <- down_csub %>% mutate(reg = "Down")

#up csub is 30, down csub is 31
csub_gene <- up_csub %>% full_join(down_csub)

write.csv(csub_gene, 'manuscript/brain/manuscript70/results/tables/csub_reorganization_MEA_genes.csv' ,row.names = F)
#########


upstatus <- sa %>% filter(logFC>0.2) 

ups <- upstatus %>% arrange(-logFC)%>% head(.,25) %>% tibble(.)
ups$symbol

down <- sa %>% filter(logFC<0.2) 

ds <- down %>% arrange(logFC)%>% head(.,25) %>% tibble(.)
ds$symbol







MEA_cort <- readRDS('manuscript/brain/manuscript70/results/RDS/limma_PostCORT_SUB_MEA.RDS') %>% 
  distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold) %>%
  filter(.,P.Value <0.05) %>% 
  left_join(., grcm38 %>% dplyr::select(symbol, entrez))%>% 
  filter(.,!is.na(entrez)) 



#################
deg.list <- list(csub,casc,sa,MEA_cort)

deg.list %>% map(~filter(., P.Value<0.05)) %>%
  map(~summarise(.,Up = sum(logFC>0.2),
                 Down = sum(logFC<0.2))) %>%
  map(~mutate(.,Total = Up + Down))


#csub
csub_cort <-MEA_cort$symbol[(MEA_cort$symbol %in% csub$symbol)] %>% as.data.frame(.)
colnames(csub_cort)<-"symbol"

csub_cort <- csub %>% filter(symbol %in% csub_cort$symbol) %>% mutate(contrast = "CON-SUB")

#casc
casc_cort <-MEA_cort$symbol[(MEA_cort$symbol %in% casc$symbol)] %>% as.data.frame(.)
colnames(casc_cort)<-"symbol"

casc_cort <- casc %>% filter(symbol %in% casc_cort$symbol)%>% mutate(contrast = "CON-ASC")

#subas
sa_cort <-MEA_cort$symbol[(MEA_cort$symbol %in% sa$symbol)] %>% as.data.frame(.)
colnames(sa_cort)<-"symbol"

sa_cort <- sa %>% filter(symbol %in% sa_cort$symbol) %>% mutate(contrast = "SUB-ASC")

#cort MEA dataframe
MEA_cort_all <- csub_cort %>% rbind(casc_cort, sa_cort)
head(MEA_cort_all)

MEA_cort_all$regStatus<- ifelse(MEA_cort_all$logFC >= 0.2, "up", "down")
table(MEA_cort_all$contrast, MEA_cort_all$regStatus)

colnames(MEA_cort_all)[2:4] <- paste0(colnames(MEA_cort_all)[2:4], "Status" )

MEA_cort_all <- MEA_cort_all %>% full_join(MEA_cort)
colnames(MEA_cort_all)[8:10] <- paste0(colnames(MEA_cort_all)[8:10], "Cort" )

MEA_cort_all$regCort<- ifelse(MEA_cort_all$logFCCort >= 0.2, "up", "down")

scort <- MEA_cort_all %>% filter(regStatus == MEA_cort_all$regCort)

table(scort$contrast, scort$regStatus)

MEA_cort_all<- anti_join(MEA_cort_all,scort)


write.csv(MEA_cort_all, 'manuscript/brain/manuscript70/results/tables/SUB_MEA_cort_DEGs.csv' ,row.names = F)
write.csv(scort, 'manuscript/brain/manuscript70/results/tables/SUB_MEA_cort_overlap_Status_DEGs_samedirection.csv' ,row.names = F)

###############
#venn Diagram 
library(vennplot)
library(systemPipeR)
deg.up <- limma_list %>% map(~filter(.,logFC >= 0.2)) %>% map(~dplyr::select(.,symbol))
up.cort <- MEA_cort %>% filter(.,logFC >= 0.2) %>% dplyr::select(.,symbol)
l.up <- list("CON-SUB" = c(t(deg.up$controlsub)),"CON-ASC" = c(t(deg.up$controlasc)),"SUB-ASC"= c(t(deg.up$subasc)), "CORT"= c(t(up.cort)))
vennset.up <- overLapper(l.up[1:4], type="vennsets")
vennPlot(vennset.up)

deg.down <- limma_list %>% map(~filter(.,logFC <= 0.2)) %>% map(~dplyr::select(.,symbol))
down.cort <- MEA_cort %>% filter(.,logFC <= 0.2) %>% dplyr::select(.,symbol)
l.down <- list("CON-SUB" = c(t(deg.down$controlsub)),"CON-ASC" = c(t(deg.down$controlasc)),"SUB-ASC"= c(t(deg.down$subasc)),"CORT"= c(t(down.cort)))
vennset.down <- overLapper(l.down[1:4], type="vennsets")
vennPlot(vennset.down)

vennPlot(list(vennset.down,vennset.up),mymain="", mysub="", colmode=2, ccol=c( "red", "blue"))

innset.down <- overLapper(l.up[1:4], type="intersects")

intersectlist(innset.down)


# mPFC data

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_mPFC_Controlsub.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

csub <- limma_list$controlsub
casc <- limma_list$controlasc
sa <- limma_list$subasc



upstatus <- sa %>% filter(logFC>0.2) 

ups <- upstatus %>% arrange(-logFC)%>% head(.,25) %>% tibble(.)
ups$symbol

down <- sa %>% filter(logFC<0.2) 

ds <- down %>% arrange(logFC)%>% head(.,25) %>% tibble(.)
ds$symbol



mPFC_cort <- readRDS('manuscript/brain/manuscript70/results/RDS/limma_PostCORT_SUB_mPFC.RDS') %>% 
  distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold) %>%
  filter(.,P.Value <0.05) %>% 
  left_join(., grcm38 %>% dplyr::select(symbol, entrez) )%>% 
  filter(.,!is.na(entrez))


#csub
csub_cort <-mPFC_cort$symbol[(mPFC_cort$symbol %in% csub$symbol)] %>% as.data.frame(.)
colnames(csub_cort)<-"symbol"

csub_cort <- csub  %>% filter(symbol %in% csub_cort$symbol) %>% mutate(contrast = "CON-SUB")

#casc
casc_cort <-mPFC_cort$symbol[(mPFC_cort$symbol %in% casc$symbol)] %>% as.data.frame(.)
colnames(casc_cort)<-"symbol"

casc_cort <- casc %>% filter(symbol %in% casc_cort$symbol)%>% mutate(contrast = "CON-ASC")

#subasc
sa_cort <-mPFC_cort$symbol[(mPFC_cort$symbol %in% sa$symbol)] %>% as.data.frame(.)
colnames(sa_cort)<-"symbol"

sa_cort <- sa %>% filter(symbol %in% sa_cort$symbol) %>% mutate(contrast = "SUB-ASC")

#cort mPFC dataframe
mPFC_cort_all <- csub_cort %>% rbind(casc_cort, sa_cort)
head(mPFC_cort_all)

mPFC_cort_all$regStatus<- ifelse(mPFC_cort_all$logFC >= 0.2, "up", "down")
table(mPFC_cort_all$contrast, mPFC_cort_all$regStatus)

colnames(mPFC_cort_all)[2:4] <- paste0(colnames(mPFC_cort_all)[2:4], "Status" )

mPFC_cort_all <- mPFC_cort_all %>% full_join(mPFC_cort)

colnames(mPFC_cort_all)[8:10] <- paste0(colnames(mPFC_cort_all)[8:10], "Cort" )

mPFC_cort_all$regCort<- ifelse(mPFC_cort_all$logFCCort >= 0.2, "up", "down")

scort <- mPFC_cort_all %>% filter(regStatus == mPFC_cort_all$regCort)

table(scort$contrast, scort$regStatus)

mPFC_cort_all<- anti_join(mPFC_cort_all,scort)
write.csv(mPFC_cort_all, 'manuscript/brain/manuscript70/results/tables/SUB_mPFC_cort_DEGs.csv' ,row.names = F)
write.csv(scort, 'manuscript/brain/manuscript70/results/tables/SUB_mPFC_cort_overlap_Status_DEGs_samedirection.csv' ,row.names = F)

###############
#venn Diagram 
library(vennplot)
library(systemPipeR)
deg.up <- limma_list %>% map(~filter(.,logFC >= 0.2)) %>% map(~dplyr::select(.,symbol))
up.cort <- mPFC_cort %>% filter(.,logFC >= 0.2) %>% dplyr::select(.,symbol)
l.up <- list("CON-SUB" = c(t(deg.up$controlsub)),"CON-ASC" = c(t(deg.up$controlasc)),"SUB-ASC"= c(t(deg.up$subasc)), "CORT"= c(t(up.cort)))
vennset.up <- overLapper(l.up[1:4], type="vennsets")
vennPlot(vennset.up)

deg.down <- limma_list %>% map(~filter(.,logFC <= 0.2)) %>% map(~dplyr::select(.,symbol))
down.cort <- mPFC_cort %>% filter(.,logFC <= 0.2) %>% dplyr::select(.,symbol)
l.down <- list("CON-SUB" = c(t(deg.down$controlsub)),"CON-ASC" = c(t(deg.down$controlasc)),"SUB-ASC"= c(t(deg.down$subasc)),"CORT"= c(t(down.cort)))
vennset.down <- overLapper(l.down[1:4], type="vennsets")
vennPlot(vennset.down)

vennPlot(list(vennset.down,vennset.up),mymain="", mysub="", colmode=2, ccol=c( "red", "blue"))


innset.down <- overLapper(l.up[1:4], type="intersects")

intersectlist(innset.down)

