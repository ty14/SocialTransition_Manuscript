library(tidyverse)


my_logFC_threshold = 0.2

#MEA data

limma_list<- readRDS("brain/results/RDS/limma_MEA_CDOM.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol,ensgene, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

cdom <- limma_list$cdom

cdes <- limma_list$cdes 

domdes <- limma_list$domdes 

#looking further into descending genes:
des <- read_csv('brain/results/tables/distruption_MEA_genes.csv')
head(des)
sd <- domdes[domdes$symbol %in% des$symbol,]


## sub 

limma_list<- readRDS("brain/results/RDS/limma_MEA_CSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

cs <- limma_list$csub

ca <- limma_list$casc 

sa <- limma_list$subasc 

#looking further into ascending genes:
asc<- read_csv('brain/results/tables/distruptionSUB_MEA_genes.csv')
head(asc)
sas <- sa[sa$symbol %in% asc$symbol,]


##social transition genes
ascx <- asc %>% select(symbol, reg_sub = reg)
desx <- des %>% select(symbol, reg_dom = reg)

dp <- ascx %>% full_join(desx)
head(dp)

dp$plasticity <- ifelse(dp$reg_dom == 'Up' & dp$reg_sub == "Up","plastic","NA")
dp$plasticity <- ifelse(dp$reg_dom == 'Down' & dp$reg_sub == "Down","plastic",dp$plasticity)
dp$plasticity1 <- ifelse(dp$reg_dom == 'Up' & is.na(dp$reg_sub),"dom",'NA')
dp$plasticity1 <- ifelse(dp$reg_dom == 'Down' & is.na(dp$reg_sub),"dom",dp$plasticity1)
dp$plasticity1 <- ifelse(dp$reg_sub == 'Up' & is.na(dp$reg_dom),"sub",dp$plasticity1)
dp$plasticity1 <- ifelse(dp$reg_sub == 'Down' & is.na(dp$reg_dom),"sub",dp$plasticity1)

#disruption
Disruption_genes <- dp %>% filter(plasticity == "plastic") %>% select(symbol)
cd <- cdes %>% select(symbol, cd_LF = logFC, cd_PV = P.Value)
ca <- ca %>% select(symbol, ca_LF = logFC, ca_PV = P.Value )

dx <- cd[cd$symbol %in% Disruption_genes$symbol,]

ax <- ca[ca$symbol %in% Disruption_genes$symbol,]

Disruption_genes2 <- dx %>% full_join(ax) 
#one gene
# symbol      cd_LF  cd_PV     ca_LF ca_PV
# 1    Mt3 -0.5105466 0.0092 -0.433208 0.047

#descenders
Dominant_genes <- dp %>% filter(plasticity1 == "dom") %>% select(symbol)
Dominant_genes2 <- cd[cd$symbol %in% Dominant_genes$symbol,]

Dominant_genes2 %>% arrange(-cd_LF) %>% head(.,10)
Dominant_genes2 %>% arrange(cd_LF) %>% head(.,10)

#ascenders
Sub_genes <- dp %>% filter(plasticity1 == "sub") %>% select(symbol)
Sub_genes2 <- ca[ca$symbol %in% Sub_genes$symbol,]

Sub_genes2 %>% arrange(-ca_LF) %>% head(.,10)
Sub_genes2 %>% arrange(ca_LF) %>% head(.,10)


#mPFC data

limma_list<- readRDS("brain/results/RDS/limma_mPFC_CDOM.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol,ensgene, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

cdom <- limma_list$cdom

cdes <- limma_list$cdes 

domdes <- limma_list$domdes 

#looking further into descending genes:
des <- read_csv('brain/results/tables/distruption_mPFC_genes.csv')
head(des)
sd <- domdes[domdes$symbol %in% des$symbol,]

## sub 

limma_list<- readRDS("brain/results/RDS/limma_mPFC_CSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

cs <- limma_list$csub

ca <- limma_list$casc 

sa <- limma_list$subasc 

asc<- read_csv('brain/results/tables/distruptionSUB_mPFC_genes.csv')
head(asc)
sas <- sa[sa$symbol %in% asc$symbol,]

##social disruption genes
ascx <- asc %>% select(symbol, reg_sub = reg)
desx <- des %>% select(symbol, reg_dom = reg)

dp <- ascx %>% full_join(desx)
head(dp)

dp$plasticity <- ifelse(dp$reg_dom == 'Up' & dp$reg_sub == "Up","plastic","NA")
dp$plasticity <- ifelse(dp$reg_dom == 'Down' & dp$reg_sub == "Down","plastic",dp$plasticity)
dp$plasticity1 <- ifelse(dp$reg_dom == 'Up' & is.na(dp$reg_sub),"dom",'NA')
dp$plasticity1 <- ifelse(dp$reg_dom == 'Down' & is.na(dp$reg_sub),"dom",dp$plasticity1)
dp$plasticity1 <- ifelse(dp$reg_sub == 'Up' & is.na(dp$reg_dom),"sub",dp$plasticity1)
dp$plasticity1 <- ifelse(dp$reg_sub == 'Down' & is.na(dp$reg_dom),"sub",dp$plasticity1)

#transitions
Disruption_genes <- dp %>% filter(plasticity == "plastic") %>% select(symbol)
cd <- cdes %>% select(symbol, cd_LF = logFC, cd_PV = P.Value)
ca <- ca %>% select(symbol, ca_LF = logFC, ca_PV = P.Value )

dx <- cd[cd$symbol %in% Disruption_genes$symbol,]

ax <- ca[ca$symbol %in% Disruption_genes$symbol,]

Disruption_genes2 <- dx %>% full_join(ax) 
#one gene
# symbol      cd_LF  cd_PV      ca_LF  ca_PV
# 1 1810055G02Rik -0.4875238 0.0118 -0.3023901 0.0358

#descenders
Dominant_genes <- dp %>% filter(plasticity1 == "dom") %>% select(symbol)
Dominant_genes2 <- cd[cd$symbol %in% Dominant_genes$symbol,]

Dominant_genes2 %>% arrange(-cd_LF) %>% head(.,10)
Dominant_genes2 %>% arrange(cd_LF) %>% head(.,10)

#ascenders
Subordinate_genes <- dp %>% filter(plasticity1 == "sub") %>% select(symbol)
Subordinate_genes2 <- ca[ca$symbol %in% Subordinate_genes$symbol,]

Subordinate_genes2 %>% arrange(-ca_LF) %>% head(.,10)
Subordinate_genes2 %>% arrange(ca_LF) %>% head(.,10)




