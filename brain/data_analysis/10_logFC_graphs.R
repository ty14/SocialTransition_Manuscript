library(limma)
library(edgeR)
library(Mus.musculus)
organism = 'org.Mm.eg.db'
library(organism, character.only = TRUE)
library(biomaRt)
library(AnnotationDbi)
library(annotables)
grcm38 # mouse genes
library(tidyverse)


my_logFC_threshold = 0.1

#MEA data

limma_list<- readRDS("brain/results/RDS/limma_MEA_CDOM.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol,ensgene, entrez))) %>% 
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


des <- up_des %>% full_join(down_des)

sd <- domdes[domdes$symbol %in% des$symbol,]

## sub 

limma_list<- readRDS("brain/results/RDS/limma_MEA_CSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

csub <- limma_list$csub
csub$reg <- ifelse(csub$logFC >= 0.2, "up", "down")

casc <- limma_list$casc 
casc$reg <- ifelse(casc$logFC >= 0.2, "up", "down")

subasc <- limma_list$subasc 
subasc$reg <- ifelse(subasc$logFC >= 0.2, "up", "down")

#looking further into ascending genes:
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

asc <- up_asc %>% full_join(down_asc)

sas <- sa[sa$symbol %in% asc$symbol,]


sd #2508
sas #1959
x <- sas[sd$symbol %in% sas$symbol,]

sdx <- sd %>% mutate(tran = "descender") %>% dplyr::select(-ensgene, -reg)
sasx <- sas %>% mutate(tran = "ascender")

tran <- sdx %>% full_join(sasx, by = 'symbol') %>% na.omit(.)

head(tran)
tail(tran)


ggplot(tran, )



head(subasc) #5056
head(domdes) #6230



dx <- domdes %>% full_join(subasc, by = "symbol") %>% na.omit(.) 

lab <- dx %>% filter(P.Value.x < 0.05 & P.Value.y < 0.05) %>% dplyr::select(symbol) 

library("ggrepel")
head(dx)

ggplot(dx, aes(logFC.x, logFC.y))+
  geom_point(color = 'grey', shape = 1, size = 2.5, alpha = 2)+ 
  geom_text_repel(data=subset(dx, P.Value.x < 0.015 & P.Value.y < 0.015),aes(logFC.x, logFC.y,label=symbol))+
  xlab("LogFC of DOM - DES") +
  ylab("LogFC of SUB - ASC") +
  theme_bw()

library(tidyverse)
cdom <- cdom %>% select(-ensgene) 

dx <- cdom %>% full_join(csub, by = "symbol") %>% na.omit(.) 

lab <- dx %>% filter(P.Value.x < 0.05 & P.Value.y < 0.05) %>% dplyr::select(symbol) 

library("ggrepel")
head(dx)

ggplot(dx, aes(logFC.x, logFC.y))+
  geom_point(color = 'grey', shape = 1, size = 2.5, alpha = 2)+ 
  geom_text_repel(data=subset(dx, P.Value.x < 0.015 & P.Value.y < 0.015),aes(logFC.x, logFC.y,label=symbol))+
  xlab("LogFC of CDOM - DOM") +
  ylab("LogFC of CSUB - SUB") +
  theme_bw()

