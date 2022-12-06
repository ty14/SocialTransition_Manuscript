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
library(goseq)
library(clusterProfiler)
library(enrichplot)
library(organism, character.only = TRUE)
library(DOSE)
library(tidyverse)
grcm38 # mouse genes


library(tidyverse)


my_logFC_threshold = 0.2

#MEA data

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_mPFC_ControlDD.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol,ensgene, entrez))) %>% 
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

up_des<- up_des %>% mutate(reg = "Up")

cddd <- cdes %>% filter(reg == "up")

dddd <- domdes %>% filter(reg == "up")

down_des <- cddd$symbol[(cddd$symbol %in% dddd$symbol)] %>% as.data.frame(.)
colnames(down_des)<-"symbol"

down_des <- down_des %>% mutate(reg = "Down")

#descenders genes
des <- up_des %>% full_join(down_des)

#mea up 270, 237
#mpfc up 75, down 74
##################################################################################################

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_mPFC_ControlSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

cs <- limma_list$controlsub
cs$reg <- ifelse(cs$logFC >= 0.2, "up", "down")

ca <- limma_list$controlasc 
ca$reg <- ifelse(ca$logFC >= 0.2, "up", "down")

sa <- limma_list$subasc 
sa$reg <- ifelse(sa$logFC >= 0.2, "up", "down")


#looking for asc upregulated genes: 
css <- ca %>% filter(reg == "down")

csa <- sa %>% filter(reg == "down")

up_asc <- css$symbol[(css$symbol %in% csa$symbol)] %>% as.data.frame(.)
colnames(up_asc)<-"symbol"

up_asc<- up_asc %>% mutate(reg = "Up")

csss <- ca %>% filter(reg == "up")

csaa <- sa %>% filter(reg == "up")

down_asc <- csss$symbol[(csss$symbol %in% csaa$symbol)] %>% as.data.frame(.)
colnames(down_asc)<-"symbol"

down_asc <- down_asc %>% mutate(reg = "Down")

#ascenders genes
asc <- up_asc %>% full_join(down_asc)

#mea up 103,  down 80
#mpfc up 57, 65 

#Social transition genes 
ascx <- asc %>% select(symbol, reg_asc = reg)
desx <- des %>% select(symbol, reg_des = reg)

x <- ascx %>% full_join(desx)
head(x)

x$transition <- ifelse(x$reg_asc == 'Up' & x$reg_des == "Up","Tran","NA")
x$transition <- ifelse(x$reg_asc == 'Down' & x$reg_des == "Down","Tran",x$transition)
x$transition1 <- ifelse(x$reg_asc == 'Up' & is.na(x$reg_des),"asc",'NA')
x$transition1 <- ifelse(x$reg_asc == 'Down' & is.na(x$reg_des),"asc",x$transition1)
x$transition1 <- ifelse(x$reg_des == 'Up' & is.na(x$reg_asc),"des",x$transition1)
x$transition1 <- ifelse(x$reg_des == 'Down' & is.na(x$reg_asc),"des",x$transition1)


Transition_genes <- x %>% filter(transition == "Tran") %>% select(symbol)

dd <- domdes %>% select(symbol, domdes_LF = logFC, domdes_PV = P.Value)
aa <- sa %>% select(symbol, subasc_LF = logFC, subasc_PV = P.Value )

dx <- dd[dd$symbol %in% Transition_genes$symbol,]

ax <- aa[aa$symbol %in% Transition_genes$symbol,]

Transition_genes2 <- dx %>% full_join(ax)


################################################################################
#Social distruption genes


#looking for cdom upregulated genes: 
cdd <- cdom %>% filter(reg == "up")

ddd <- cdes %>% filter(reg == "up")

up_dom <- cdd$symbol[(cdd$symbol %in% ddd$symbol)] %>% as.data.frame(.)
colnames(up_dom)<-"symbol"

up_dom<- up_dom %>% mutate(reg = "Up") #115, 83 (mpfc)

cddd <- cdom %>% filter(reg == "down")

dddd <- cdes %>% filter(reg == "down")

down_dom <- cddd$symbol[(cddd$symbol %in% dddd$symbol)] %>% as.data.frame(.)
colnames(down_dom)<-"symbol"

down_dom <- down_dom %>% mutate(reg = "Down") #98,67 (mpfc)

## distruption  cdom genes
dom <- up_dom %>% full_join(down_dom)

####################
#looking for cdom upregulated genes: 
css <- cs %>% filter(reg == "up")

csa <- ca %>% filter(reg == "up")

up_sub <- css$symbol[(css$symbol %in% csa$symbol)] %>% as.data.frame(.)
colnames(up_sub)<-"symbol"

up_sub<- up_sub %>% mutate(reg = "Up") #30, 113 (mpfc)

csss <- cs %>% filter(reg == "down")

csaa <- ca %>% filter(reg == "down")

down_sub <- csss$symbol[(csss$symbol %in% csaa$symbol)] %>% as.data.frame(.)
colnames(down_sub)<-"symbol"

down_sub <- down_sub %>% mutate(reg = "Down") #31, 111 (mpfc)

## distruption  csub genes
sub<- up_sub %>% full_join(down_sub)


domx <- dom %>% select(symbol, reg_dom = reg)
subx <- sub %>% select(symbol, reg_sub = reg)

dp <- domx %>% full_join(subx)
head(dp)

dp$plasticity <- ifelse(dp$reg_dom == 'Up' & dp$reg_sub == "Up","plastic","NA")
dp$plasticity <- ifelse(dp$reg_dom == 'Down' & dp$reg_sub == "Down","plastic",dp$plasticity)
dp$plasticity1 <- ifelse(dp$reg_dom == 'Up' & is.na(dp$reg_sub),"dom",'NA')
dp$plasticity1 <- ifelse(dp$reg_dom == 'Down' & is.na(dp$reg_sub),"dom",dp$plasticity1)
dp$plasticity1 <- ifelse(dp$reg_sub == 'Up' & is.na(dp$reg_dom),"sub",dp$plasticity1)
dp$plasticity1 <- ifelse(dp$reg_sub == 'Down' & is.na(dp$reg_dom),"sub",dp$plasticity1)

plas <- dp %>% filter(plasticity == "plastic") %>% select(symbol)
plas 
#Mt3 mea
#mpfc: 
#Entpd2
# Dph7
# 1810055G02Rik

dom <- dp %>% filter(plasticity1 == "dom") %>% select(symbol)
dom 
#210 mea
#147 mpfc

sub <- dp %>% filter(plasticity1 == "sub") %>% select(symbol)
sub 
#58 mea
#221 mpfc 


cdom <- cdom %>% select(symbol, cdom_LF = logFC, cdom_PV = P.Value)
csub <- cs %>% select(symbol, csub_LF = logFC, csub_PV = P.Value )


cdom[cdom$symbol %in% plas$symbol,] 
##mea Mt3 up regulated in doms
## mpfc: dph7/1810055G02Rik up in dom, entpd2 up in control

csub[csub$symbol %in% plas$symbol,]
##mea Mt3 up regulated in subs
## mpfc: dph7/1810055G02Rik up in sub, entpd2 up in control



