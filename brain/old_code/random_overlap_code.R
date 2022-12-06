library(tidyverse)


my_logFC_threshold = 0.2

#MEA data

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_MEA_ControlDD.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
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
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

cs <- limma_list$controlsub

ca <- limma_list$controlasc 

sa <- limma_list$subasc 


asc <- read_csv('manuscript/brain/manuscript70/results/tables/ASC_transition_mPFC_genes.csv')
head(asc)


des <- read_csv('manuscript/brain/manuscript70/results/tables/descenders_tranisiton_mPFC_genes.csv')
head(des)


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



head(x)

Descender_genes <- x %>% filter(transition1 == "des") %>% select(symbol)
Descender_genes2 <- dd[dd$symbol %in% Descender_genes$symbol,]

Descender_genes2 %>% filter(domdes_LF > -0.2)

Ascension_genes <- x %>% filter(transition1 == "asc") %>% select(symbol)
Ascension_genes2 <- aa[aa$symbol %in% Ascension_genes$symbol,]



write.csv(Ascension_genes2, "manuscript/brain/manuscript70/results/tables/ASC_MEA_genes.csv",row.names = F)

write.csv(Descender_genes2, "manuscript/brain/manuscript70/results/tables/DES_MEA_genes.csv",row.names = F)

write.csv(Transition_genes2, "manuscript/brain/manuscript70/results/tables/transition_MEA_genes.csv",row.names = F)


########
# plasticity

dom <- read_csv('manuscript/brain/manuscript70/results/tables/dominant_distruption_MEA_genes.csv')
head(dom)

sub <- read_csv('manuscript/brain/manuscript70/results/tables/subordinate_distruption_MEA_genes.csv')
head(sub)

d <- sub %>% filter(reg == "Down")
p <- sub %>% filter(reg == "Up")


domx <- dom %>% select(symbol, reg_dom = reg)
subx <- sub %>% select(symbol, reg_sub = reg)

dp <- domx %>% full_join(subx)
head(dp)


dp$plasticity <- ifelse(dp$reg_dom == 'Up' & dp$reg_sub == "Up","plastic","NA")
dp$plasticity <- ifelse(dp$reg_dom == 'Down' & dp$reg_sub == "Down","plastic",dp$plasticity)


# dp$plasticity1 <- ifelse(dp$dom_reg == 'Up' & is.na(dp$sub_reg),"dom",'NA')
# dp$plasticity1 <- ifelse(dp$dom_reg == 'Down' & is.na(dp$sub_reg),"dom",dp$plasticity1)
# dp$plasticity1 <- ifelse(dp$sub_reg == 'Up' & is.na(dp$dom_reg),"sub",dp$plasticity1)
# dp$plasticity1 <- ifelse(dp$sub_reg == 'Down' & is.na(dp$dom_reg),"sub",dp$plasticity1)

dp$plasticity1 <- ifelse(dp$reg_sub == 'Down' & dp$reg_dom == "Up","dom","NA")
dp$plasticity1 <- ifelse(dp$reg_sub == 'Up' & dp$reg_dom == "Down","sub",dp$plasticity1)
table(dp$plasticity1)

dp %>% filter(plasticity == 'plastic')
dp %>% filter(plasticity1 == 'dom')
dp %>% filter(plasticity1 == 'sub')


dp$reg_sub ='Down' %in% dp$reg_dom = "Up"
