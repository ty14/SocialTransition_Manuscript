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
des <- read_csv('brain/results/tables/descender_MEA_genes.csv')
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
asc<- read_csv('brain/results/tables/ascender_MEA_genes.csv')
head(asc)
sas <- sa[sa$symbol %in% asc$symbol,]


##social transition genes
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

#transitions
Transition_genes <- x %>% filter(transition == "Tran") %>% select(symbol)
dd <- domdes %>% select(symbol, domdes_LF = logFC, domdes_PV = P.Value)
aa <- sa %>% select(symbol, subasc_LF = logFC, subasc_PV = P.Value )

dx <- dd[dd$symbol %in% Transition_genes$symbol,]

ax <- aa[aa$symbol %in% Transition_genes$symbol,]

Transition_genes2 <- dx %>% full_join(ax) #5 up , 3 down  

#descenders
Descenders_genes <- x %>% filter(transition1 == "des") %>% select(symbol)
Descenders_genes2 <- dd[dd$symbol %in% Descenders_genes$symbol,]

Descenders_genes2 %>% arrange(-domdes_LF) %>% head(.,10)
Descenders_genes2 %>% arrange(domdes_LF) %>% head(.,10)

#ascenders
Ascenders_genes <- x %>% filter(transition1 == "asc") %>% select(symbol)
Ascenders_genes2 <- aa[aa$symbol %in% Ascenders_genes$symbol,]

Ascenders_genes2 %>% arrange(-subasc_LF) %>% head(.,10)
Ascenders_genes2 %>% arrange(subasc_LF) %>% head(.,10)


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
des <- read_csv('brain/results/tables/descender_mPFC_genes.csv')
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

asc<- read_csv('brain/results/tables/ascender_mPFC_genes.csv')
head(asc)
sas <- sa[sa$symbol %in% asc$symbol,]

##social transition genes
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

#transitions
Transition_genes <- x %>% filter(transition == "Tran") %>% select(symbol)
dd <- domdes %>% select(symbol, domdes_LF = logFC, domdes_PV = P.Value)
aa <- sa %>% select(symbol, subasc_LF = logFC, subasc_PV = P.Value )

dx <- dd[dd$symbol %in% Transition_genes$symbol,]

ax <- aa[aa$symbol %in% Transition_genes$symbol,]

Transition_genes2 <- dx %>% full_join(ax) #no genes

#descenders
Descenders_genes <- x %>% filter(transition1 == "des") %>% select(symbol)
Descenders_genes2 <- dd[dd$symbol %in% Descenders_genes$symbol,]

Descenders_genes2 %>% arrange(-domdes_LF) %>% head(.,10)
Descenders_genes2 %>% arrange(domdes_LF) %>% head(.,10)

#ascenders
Ascenders_genes <- x %>% filter(transition1 == "asc") %>% select(symbol)
Ascenders_genes2 <- aa[aa$symbol %in% Ascenders_genes$symbol,]

Ascenders_genes2 %>% arrange(-subasc_LF) %>% head(.,10)
Ascenders_genes2 %>% arrange(subasc_LF) %>% head(.,10)



