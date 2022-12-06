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

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_MEA_ControlDD.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

cdom <- limma_list$controldom
  cdom$reg <- ifelse(cdom$logFC >= 0.2, "up", "down")

cdes <- limma_list$controldes 
  cdes$reg <- ifelse(cdes$logFC >= 0.2, "up", "down")


domdes <- limma_list$domdes 
  domdes$reg <- ifelse(domdes$logFC >= 0.2, "up", "down")


#looking for des upregulated genes: 
cdd <- cdes %>% filter(reg == "down")

ddd <- domdes %>% filter(reg == "down")

up_des <- cdd$symbol[(cdd$symbol %in% ddd$symbol)] %>% as.data.frame(.)
colnames(up_des)<-"symbol"

up_des <- up_des %>% mutate(reg = "Up")

cddd <- cdes %>% filter(reg == "up")

dddd <- domdes %>% filter(reg == "up")

down_des <- cddd$symbol[(cddd$symbol %in% dddd$symbol)] %>% as.data.frame(.)
colnames(down_des)<-"symbol"

down_des <- down_des %>% mutate(reg = "Down")

#up descenders is 270, down decenders is 237
des <- up_des %>% full_join(down_des)


write.csv(des, 'manuscript/brain/manuscript70/results/tables/descenders_tranisiton_MEA_genes.csv' ,row.names = F)

## looking for genes different from CDOM: 
cdd <- cdes %>% filter(reg == "up")

cddom <- cdom %>% filter(reg == "up")

up_cdom <- cdd$symbol[(cdd$symbol %in% cddom$symbol)] %>% as.data.frame(.)
colnames(up_cdom)<-"symbol"

up_cdom <- up_cdom %>% mutate(reg = "Up")

cddd <- cdes %>% filter(reg == "down")

cdddom <- cdom %>% filter(reg == "down")

down_cdom <- cddd$symbol[(cddd$symbol %in% cdddom$symbol)] %>% as.data.frame(.)
colnames(down_cdom)<-"symbol"

down_cdom <- down_cdom %>% mutate(reg = "Down")

#up cdom is 115, down cdom is 98
cdom_gene <- up_cdom %>% full_join(down_cdom)

write.csv(cdom_gene, 'manuscript/brain/manuscript70/results/tables/cdom_reorganization_MEA_genes.csv' ,row.names = F)

#########

upstatus <- domdes %>% filter(logFC>0.2) 

ups <- upstatus %>% arrange(-logFC)%>% head(.,25) %>% tibble(.)
ups$symbol

down <- domdes %>% filter(logFC<0.2) 

ds <- down %>% arrange(logFC)%>% head(.,25) %>% tibble(.)
ds$symbol


MEA_cort <- readRDS('manuscript/brain/manuscript70/results/RDS/limma_PostCORT_DOM_MEA.RDS') %>% 
  distinct(.) %>% 
 filter(.,abs(logFC) >= my_logFC_threshold) %>%
 filter(.,P.Value <0.05) %>% 
left_join(., grcm38 %>% dplyr::select(symbol, entrez))%>% 
  filter(.,!is.na(entrez)) 


#################
deg.list %>% map(~filter(., P.Value<0.05)) %>%
  map(~summarise(.,Up = sum(logFC>0.2),
                 Down = sum(logFC<0.2))) %>%
  map(~mutate(.,Total = Up + Down))


#cdom
cdom_cort <-MEA_cort$symbol[(MEA_cort$symbol %in% cdom$symbol)] %>% as.data.frame(.)
colnames(cdom_cort)<-"symbol"

cdom_cort <- cdom %>% filter(symbol %in% cdom_cort$symbol) %>% mutate(contrast = "CON-DOM")


#cdes
cdes_cort <-MEA_cort$symbol[(MEA_cort$symbol %in% cdes$symbol)] %>% as.data.frame(.)
colnames(cdes_cort)<-"symbol"

cdes_cort <- cdes %>% filter(symbol %in% cdes_cort$symbol)%>% mutate(contrast = "CON-DES")

#domdes
dd_cort <-MEA_cort$symbol[(MEA_cort$symbol %in% domdes$symbol)] %>% as.data.frame(.)
colnames(dd_cort)<-"symbol"

dd_cort <- domdes %>% filter(symbol %in% dd_cort$symbol) %>% mutate(contrast = "DOM-DES")

#cort MEA dataframe
MEA_cort_all <- cdom_cort %>% rbind(cdes_cort, dd_cort)
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

write.csv(MEA_cort_all, 'manuscript/brain/manuscript70/results/tables/DOM_MEA_cort_DEGs.csv' ,row.names = F)
write.csv(scort, 'manuscript/brain/manuscript70/results/tables/DOM_MEA_cort_overlap_Status_DEGs_samedirection.csv' ,row.names = F)

###############
#venn Diagram 
library(vennplot)
library(systemPipeR)
deg.up <- limma_list %>% map(~filter(.,logFC >= 0.2)) %>% map(~dplyr::select(.,symbol))
up.cort <- MEA_cort %>% filter(.,logFC >= 0.2) %>% dplyr::select(.,symbol)
l.up <- list("CON-DOM" = c(t(deg.up$controldom)),"CON-DES" = c(t(deg.up$controldes)),"DOM-DES"= c(t(deg.up$domdes)), "CORT"= c(t(up.cort)))
vennset.up <- overLapper(l.up[1:4], type="vennsets")
vennPlot(vennset.up)

deg.down <- limma_list %>% map(~filter(.,logFC <= 0.2)) %>% map(~dplyr::select(.,symbol))
down.cort <- MEA_cort %>% filter(.,logFC <= 0.2) %>% dplyr::select(.,symbol)
l.down <- list("CON-DOM" = c(t(deg.down$controldom)),"CON-DES" = c(t(deg.down$controldes)),"DOM-DES"= c(t(deg.down$domdes)),"CORT"= c(t(down.cort)))
vennset.down <- overLapper(l.down[1:4], type="vennsets")
vennPlot(vennset.down)

vennPlot(list(vennset.down,vennset.up),mymain="", mysub="", colmode=2, ccol=c( "red", "blue"))



# mPFC data

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_mPFC_ControlDD.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

cdom <- limma_list$controldom
cdes <- limma_list$controldes 
domdes <- limma_list$domdes 



upstatus <- cdes %>% filter(logFC>0.2) 

ups <- upstatus %>% arrange(-logFC)%>% head(.,25) %>% tibble(.)
ups$symbol

down <- cdes %>% filter(logFC<0.2) 

ds <- down %>% arrange(logFC)%>% head(.,25) %>% tibble(.)
ds$symbol


mPFC_cort <- readRDS('manuscript/brain/manuscript70/results/RDS/limma_PostCORT_DOM_mPFC.RDS') %>% 
  distinct(.) %>% 
  filter(.,abs(logFC) >= my_logFC_threshold) %>%
  filter(.,P.Value <0.05) %>% 
  left_join(., grcm38 %>% dplyr::select(symbol, entrez) )%>% 
  filter(.,!is.na(entrez))


#cdom
cdom_cort <-mPFC_cort$symbol[(mPFC_cort$symbol %in% cdom$symbol)] %>% as.data.frame(.)
colnames(cdom_cort)<-"symbol"

cdom_cort <- cdom  %>% filter(symbol %in% cdom_cort$symbol) %>% mutate(contrast = "CON-DOM")


#cdes
cdes_cort <-mPFC_cort$symbol[(mPFC_cort$symbol %in% cdes$symbol)] %>% as.data.frame(.)
colnames(cdes_cort)<-"symbol"

cdes_cort <- cdes %>% filter(symbol %in% cdes_cort$symbol)%>% mutate(contrast = "CON-DES")

#domdes
dd_cort <-mPFC_cort$symbol[(mPFC_cort$symbol %in% domdes$symbol)] %>% as.data.frame(.)
colnames(dd_cort)<-"symbol"

dd_cort <- domdes %>% filter(symbol %in% dd_cort$symbol) %>% mutate(contrast = "DOM-DES")

#cort mPFC dataframe
mPFC_cort_all <- cdom_cort %>% rbind(cdes_cort, dd_cort)
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
write.csv(mPFC_cort_all, 'manuscript/brain/manuscript70/results/tables/DOM_mPFC_cort_DEGs.csv' ,row.names = F)
write.csv(scort, 'manuscript/brain/manuscript70/results/tables/DOM_mPFC_cort_overlap_Status_DEGs_samedirection.csv' ,row.names = F)

###############
#venn Diagram 
library(vennplot)
library(systemPipeR)
deg.up <- limma_list %>% map(~filter(.,logFC >= 0.2)) %>% map(~dplyr::select(.,symbol))
up.cort <- mPFC_cort %>% filter(.,logFC >= 0.2) %>% dplyr::select(.,symbol)
l.up <- list("CON-DOM" = c(t(deg.up$controldom)),"CON-DES" = c(t(deg.up$controldes)),"DOM-DES"= c(t(deg.up$domdes)), "CORT"= c(t(up.cort)))
vennset.up <- overLapper(l.up[1:4], type="vennsets")
vennPlot(vennset.up)

deg.down <- limma_list %>% map(~filter(.,logFC <= 0.2)) %>% map(~dplyr::select(.,symbol))
down.cort <- mPFC_cort %>% filter(.,logFC <= 0.2) %>% dplyr::select(.,symbol)
l.down <- list("CON-DOM" = c(t(deg.down$controldom)),"CON-DES" = c(t(deg.down$controldes)),"DOM-DES"= c(t(deg.down$domdes)),"CORT"= c(t(down.cort)))
vennset.down <- overLapper(l.down[1:4], type="vennsets")
vennPlot(vennset.down)

vennPlot(list(vennset.down,vennset.up),mymain="", mysub="", colmode=2, ccol=c( "red", "blue"))


innset.down <- overLapper(l.down[1:4], type="intersects")

intersectlist(innset.down)

