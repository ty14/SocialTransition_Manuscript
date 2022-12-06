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


my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_MEA_ControlDD.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez, description))) %>% 
  map(~filter(.,!is.na(entrez))) 

cdom <- limma_list$controldom
cdes <- limma_list$controldes 
domdes <- limma_list$domdes 


cdom %>% filter(grepl("Chr", symbol))
cdes %>% filter(grepl("Chr", symbol))

domdes %>% filter(grepl("Chr", symbol))

#control dominant
upcdom <- cdom %>% 
  filter(logFC>0.2) 

upcd <- upcdom %>% arrange(-logFC) %>% head(25) %>% tibble(.)

downcdom <- cdom %>% 
  filter(logFC<0.2) 

downcd <- downcdom %>% arrange(logFC) %>% head(25) %>% tibble(.)

#control descender
upcdes <- cdes %>% 
  filter(logFC>0.2) 
upde <- upcdes %>% arrange(-logFC) %>% head(25) %>% tibble(.)

downcdes <- cdes %>% 
  filter(logFC<0.2) 

downde <- downcdes %>% arrange(logFC) %>% head(25) %>%  tibble(.)

#dominant descender
updomdes <- domdes %>% 
  filter(logFC>0.2) 
updd <- updomdes %>% arrange(-logFC) %>% head(25) %>%  tibble(.)

downdomdes <- domdes %>% 
  filter(logFC<0.2) 

ddd <- downdomdes %>% arrange(logFC) %>% head(25) %>%  tibble(.)



ex <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_MeA_ControlDD.RDS")
head(ex)
x <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

c1 <- x %>% filter(symbol == "Chrm1")
c2 <- x %>% filter(symbol == "Chrm2")
c3 <- x %>% filter(symbol == "Chrm3")
ch2 <- x %>% filter(symbol == "Chrnb2")
ch4 <- x %>% filter(symbol == "Chrna4")
ch6 <- x %>% filter(symbol == "Chrna6")
l <- x %>% filter(symbol == "Lhx8")
g <- x %>% filter(symbol == "Gbx2")
s18 <- x %>% filter(symbol == "Slc18a3")
s10 <- x %>% filter(symbol == "Slc10a4")


gh <- x %>% filter(symbol == "Ghrh")
th <- x %>% filter(symbol == "Trh")
ch <- x %>% filter(grepl("crh", symbol))

ht2 <-c1 %>% rbind(c2, c3,ch2,ch4,ch6,l,g,s18,s10,gh,th,ch)

id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

colnames(ht2)
ht2 <- ht2 %>% pivot_longer(cols = 2:21, names_to = "ids")

p <- ht2 %>% full_join(id)

p$group <- factor(p$group, levels = c("CDOM", "DOM", "DES"))

p <- p %>% filter(symbol != "Uqcrh")

p$symbol<- factor(p$symbol, levels = c("Chrm1", "Chrm2", "Chrm3","Chrnb2", "Chrna4","Chrna6", "Lhx8", "Gbx2", 
                                      "Slc10a4", "Slc18a3", "Ghrh", "Trh"))

source('functions/geom_boxjitter.R')

p1 <- ggplot(p, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  facet_wrap(~symbol,scales="free_y")+
  ylab("70 min expression") +
  theme_bw()+
  theme(legend.position = "none")
p1


ggsave("manuscript/brain/manuscript70/results/results_figures/DOM70_ACHgenes.png", p1, 
       width = 14, height = 10,dpi = 150)

ex <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_MeA_ControlSUB.RDS")
head(ex)


x <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

c1 <- x %>% filter(symbol == "Chrm1")
c2 <- x %>% filter(symbol == "Chrm2")
c3 <- x %>% filter(symbol == "Chrm3")
ch2 <- x %>% filter(symbol == "Chrnb2")
ch4 <- x %>% filter(symbol == "Chrna4")
ch6 <- x %>% filter(symbol == "Chrna6")
l <- x %>% filter(symbol == "Lhx8")
g <- x %>% filter(symbol == "Gbx2")
s18 <- x %>% filter(symbol == "Slc18a3")
s10 <- x %>% filter(symbol == "Slc10a4")


gh <- x %>% filter(symbol == "Ghrh")
th <- x %>% filter(symbol == "Trh")
ch <- x %>% filter(grepl("crh", symbol))




ht2 <-c1 %>% rbind(c2, c3,ch2,ch4,ch6,l,g,s18,s10,gh,th,ch)

id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

colnames(ht2)
ht2 <- ht2 %>% pivot_longer(cols = 2:21, names_to = "ids")

p <- ht2 %>% full_join(id)

p$group <- factor(p$group, levels = c("CSUB", "SUB", "ASC"))

p <- p %>% filter(symbol != "Uqcrh")

p$symbol<- factor(p$symbol, levels = c("Chrm1", "Chrm2", "Chrm3","Chrnb2", "Chrna4","Chrna6", "Lhx8", "Gbx2", 
                                       "Slc10a4", "Slc18a3", "Ghrh", "Trh"))

source('functions/geom_boxjitter.R')

p1 <- ggplot(p, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  facet_wrap(~symbol,scales="free_y")+
  ylab("70 min expression") +
  theme_bw()+
  theme(legend.position = "none")
p1


ggsave("manuscript/brain/manuscript70/results/results_figures/SUB70_ACHgenes.png", p1, 
       width = 14, height = 10,dpi = 150)

# interesting genes


my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_MEA_ControlDD.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez, description))) %>% 
  map(~filter(.,!is.na(entrez))) 


y1a <- limma_list$controldom %>% mutate(contrast = "control-dom")

y2a <- limma_list$controldes %>% mutate(contrast = "control-des")

y3a <- limma_list$domdes %>% mutate(contrast = "dom - des")

# mPFC_dom_DEGS <- y1a %>% rbind(y2a, y3a)

 MEA_dom_DEGs <- y1a %>% rbind(y2a, y3a)



limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_MEA_ControlSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez, description))) %>% 
  map(~filter(.,!is.na(entrez))) 


y1b <- limma_list$controlsub %>% mutate(contrast = "control-sub")

y2b <- limma_list$controlasc %>% mutate(contrast = "control-asc")

y3b <- limma_list$subasc %>% mutate(contrast = "sub - asc")

 # mPFC_sub_DEGs <- y1b %>% rbind(y2b, y3b)

  MEA_sub_DEGs <- y1b %>% rbind(y2b, y3b)


ex <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_MeA_ControlDD.RDS")
# head(ex)


dom_mea_exp <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

dom_m_id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

# sub_m_id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

ex <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_mPFC_ControlDD.RDS")
# head(ex)


dom_mpfc_exp <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

dom_p_id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)
 # sub_p_id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

#Hkdc1
mPFC_dom_DEGS %>% filter(symbol== "Hkdc1") # higher in doms and des
mPFC_sub_DEGs %>% filter(symbol== "Hkdc1") #higher in asc


#Hcrt 
MEA_dom_DEGs %>% filter(symbol== "Hcrt") # down in descenders 
MEA_sub_DEGs %>% filter(symbol== "Hcrt") #none

# Ghrh 
MEA_dom_DEGs %>% filter(symbol== "Ghrh") # down in descenders 
MEA_sub_DEGs %>% filter(symbol== "Ghrh") #none


# Igf1
MEA_dom_DEGs %>% filter(symbol == "Igf1")# very high in dom
MEA_sub_DEGs %>% filter(symbol== 'Igf1') # none
mPFC_dom_DEGS %>% filter(symbol == "Igf1")# none
mPFC_sub_DEGs %>% filter(symbol== 'Igf1') #none 


# Trh 
MEA_dom_DEGs %>% filter(symbol== "Trh") # down in doms
MEA_sub_DEGs %>% filter(symbol== "Trh") #none


# Chrm1
MEA_dom_DEGs %>% filter(symbol== "Chrm1") # up in des 
MEA_sub_DEGs %>% filter(symbol== "Chrm1") # up in asc 

# Chrm2
MEA_dom_DEGs %>% filter(symbol== "Chrm2") # up in doms 
MEA_sub_DEGs %>% filter(symbol== "Chrm2") # up in ascenders (aggression?)

# Chrm3
MEA_dom_DEGs %>% filter(symbol== "Chrm3") # none
MEA_sub_DEGs %>% filter(symbol== "Chrm3") # high in sub

# Chrna4
MEA_dom_DEGs %>% filter(symbol== "Chrna4") #lower in doms
MEA_sub_DEGs %>% filter(symbol== "Chrna4") # none 


# Chrna6
MEA_dom_DEGs %>% filter(symbol== "Chrna6")
MEA_sub_DEGs %>% filter(symbol== "Chrna6") # none 

# Chat 
MEA_dom_DEGs %>% filter(symbol== "Chat") #just dom - des up in dom 
MEA_sub_DEGs %>% filter(symbol== "Chat") # none 

# Slc18a3
MEA_dom_DEGs %>% filter(symbol== "Slc18a3") #very high in doms
MEA_sub_DEGs %>% filter(symbol== "Slc18a3") # none 

# Slc10a4
MEA_dom_DEGs %>% filter(symbol== "Slc10a4") #very high in doms
MEA_sub_DEGs %>% filter(symbol== "Slc10a4") # none 

# Lhx8 
MEA_dom_DEGs %>% filter(symbol== "Lhx8")#very high in doms
MEA_sub_DEGs %>% filter(symbol== "Lhx8") # none 

# Lhx1
MEA_dom_DEGs %>% filter(symbol== "Lhx1") # none
MEA_sub_DEGs %>% filter(symbol== "Lhx1") # lower in reorg subs


# Crym 
MEA_dom_DEGs %>% filter(symbol== "Crym") # very high in doms
mPFC_dom_DEGS %>% filter(symbol== "Crym") # none 
MEA_sub_DEGs %>% filter(symbol== "Crym") # high in asc 
mPFC_sub_DEGs %>% filter(symbol== "Crym") # high in asc  both in MEA and mPFC

# Mybpc1 
MEA_dom_DEGs %>% filter(symbol== "Mybpc1")# very high in doms
mPFC_dom_DEGS %>% filter(symbol== "Mybpc1")#none
MEA_sub_DEGs %>% filter(symbol== "Mybpc1") # none 
mPFC_sub_DEGs %>% filter(symbol== "Mybpc1") # high in ascenders 

# Mas1 
MEA_dom_DEGs %>% filter(symbol== "Mas1")# very high in doms
mPFC_dom_DEGS %>% filter(symbol== "Mas1")# higher in des 

MEA_sub_DEGs %>% filter(symbol== "Mas1") # none 
mPFC_sub_DEGs %>% filter(symbol== "Mas1") # none


# Foxo6
MEA_dom_DEGs %>% filter(symbol== "Foxo6")#  high in des
MEA_sub_DEGs %>% filter(symbol== "Foxo6") # very high in asc

# Cbln2
MEA_dom_DEGs %>% filter(symbol== "Cbln2")# high in doms
MEA_sub_DEGs %>% filter(symbol== "Cbln2") # very hihg in subs
mPFC_dom_DEGS %>% filter(symbol== "Cbln2")# high in doms
mPFC_sub_DEGs %>% filter(symbol== "Cbln2") # very hihg in subs



# Ogn
MEA_dom_DEGs %>% filter(symbol== "Ogn")# none
MEA_sub_DEGs %>% filter(symbol== "Ogn") # lower in subs  

# Serpina3i
mPFC_dom_DEGS %>% filter(symbol== "Serpina3i")# high in control
mPFC_sub_DEGs %>% filter(symbol== "Serpina3i") # very high in control 

# Serpina1c
MEA_dom_DEGs %>% filter(symbol== "Serpina1c")# none
MEA_sub_DEGs %>% filter(symbol== "Serpina1c") # lower in subs 

# Sim1
MEA_dom_DEGs %>% filter(symbol== "Sim1")# none
MEA_sub_DEGs %>% filter(symbol== "Sim1") # high in sub

# Rgs13 defeat stress 
MEA_dom_DEGs %>% filter(symbol== "Rgs13")# 
mPFC_dom_DEGS %>% filter(symbol== "Rgs13")# higher in des 
MEA_sub_DEGs %>% filter(symbol== "Rgs13") # very high in subs
mPFC_sub_DEGs %>% filter(symbol== "Rgs13") # high in controls

# Tmem40
MEA_dom_DEGs %>% filter(symbol== "Tmem40")# none
MEA_sub_DEGs %>% filter(symbol== "Tmem40") # very high in subs

# Ggt1
MEA_dom_DEGs %>% filter(symbol== "Ggt1")# none
MEA_sub_DEGs %>% filter(symbol== "Ggt1") # very high in subs

# Cyp26b1 maintaining status
MEA_dom_DEGs %>% filter(symbol== "Cyp26b1")# high in dom
MEA_sub_DEGs %>% filter(symbol== "Cyp26b1") # very high in subs

# Slc6a3 dopamine
MEA_dom_DEGs %>% filter(symbol== "Slc6a3")# very high in doms
MEA_sub_DEGs %>% filter(symbol== "Slc6a3") # very high in subs

# Neurod6
MEA_dom_DEGs %>% filter(symbol== "Neurod6")# none
mPFC_dom_DEGS %>% filter(symbol== "Neurod6")# high in dom
MEA_sub_DEGs %>% filter(symbol== "Neurod6") # very high in subs
mPFC_sub_DEGs %>% filter(symbol== "Neurod6") # none

# Esrrb
MEA_dom_DEGs %>% filter(symbol== "Esrrb")# very high in doms
MEA_sub_DEGs %>% filter(symbol== "Esrrb") # high in subs


# Esr1
MEA_dom_DEGs %>% filter(symbol == 'Esr1')# high in des
MEA_sub_DEGs %>% filter(symbol== 'Esr1') # none

# Ar
MEA_dom_DEGs %>% filter(symbol == "Ar")# high in control
MEA_sub_DEGs %>% filter(symbol== 'Ar') # none


# Pgrmc1
MEA_dom_DEGs %>% filter(symbol == "Pgrmc1")# high in des
MEA_sub_DEGs %>% filter(symbol== 'Pgrmc1') # none
mPFC_dom_DEGS %>% filter(symbol == "Pgrmc1")# none
mPFC_sub_DEGs %>% filter(symbol== 'Pgrmc1') # high in asc


# Pnmt
MEA_dom_DEGs %>% filter(symbol== "Pnmt")# none
MEA_sub_DEGs %>% filter(symbol== "Pnmt") # high in subs


# Arhgap36 
MEA_dom_DEGs %>% filter(symbol== "Arhgap36")# high in control dom
mPFC_dom_DEGS %>% filter(symbol== "Arhgap36")# none
MEA_sub_DEGs %>% filter(symbol== "Arhgap36") # none
mPFC_sub_DEGs %>% filter(symbol== "Arhgap36") # high in controls and subs

# Nxph4
MEA_dom_DEGs %>% filter(symbol== "Nxph4")# high in control
mPFC_dom_DEGS %>% filter(symbol== "Nxph4")# none
MEA_sub_DEGs %>% filter(symbol== "Nxph4") # high in subs
mPFC_sub_DEGs %>% filter(symbol== "Nxph4") # high in controls

# Tmem232
mPFC_dom_DEGS %>% filter(symbol== "Tmem232")# none
mPFC_sub_DEGs %>% filter(symbol== "Tmem232") # very high in reorg subs and asc

# Ttr
MEA_dom_DEGs %>% filter(symbol== "Ttr")# high in dom
mPFC_dom_DEGS %>% filter(symbol== "Ttr")# none 
MEA_sub_DEGs %>% filter(symbol== "Ttr") # none
mPFC_sub_DEGs %>% filter(symbol== "Ttr") # high in subs and asc

# Oxt
mPFC_dom_DEGS %>% filter(symbol== "Oxt")# high in control
MEA_dom_DEGs %>% filter(symbol== "Oxt")
MEA_sub_DEGs %>% filter(symbol== "Oxt")# high in control
mPFC_sub_DEGs %>% filter(symbol== "Oxt")

# Dmrtb1
MEA_dom_DEGs %>% filter(symbol== "Dmrtb1")# high in control and dom
MEA_sub_DEGs %>% filter(symbol== "Dmrtb1") # none

# Dmkn
mPFC_dom_DEGS %>% filter(symbol== "Dmkn")# high in control and dom 
mPFC_sub_DEGs %>% filter(symbol== "Dmkn") # none

# Mpeg1
MEA_dom_DEGs %>% filter(symbol== "Mpeg1")# high in dom
mPFC_sub_DEGs %>% filter(symbol== "Mpeg1") # high in control and sub


# Il1b
mPFC_dom_DEGS %>% filter(symbol== "Il1b")# high in des
mPFC_sub_DEGs %>% filter(symbol== "Il1b") #none


# Mybpc1
MEA_dom_DEGs %>% filter(symbol== "Mybpc1")# low in control dom
mPFC_dom_DEGS %>% filter(symbol== "Mybpc1")# none
MEA_sub_DEGs %>% filter(symbol== "Mybpc1") # none
mPFC_sub_DEGs %>% filter(symbol== "Mybpc1") # low in subs

# cd14
mPFC_dom_DEGS %>% filter(symbol== "Cd14")# high in dom and des
mPFC_sub_DEGs %>% filter(symbol== "Cd14") # none

# cdkl1
MEA_dom_DEGs %>% filter(symbol== "Cdkl1")# none
MEA_sub_DEGs %>% filter(symbol== "Cdkl1") # high cont


# Ndst4
MEA_dom_DEGs %>% filter(symbol== "Ndst4")# none
mPFC_dom_DEGS %>% filter(symbol== "Ndst4")# high in dom
MEA_sub_DEGs %>% filter(symbol== "Ndst4") # high in control
mPFC_sub_DEGs %>% filter(symbol== "Ndst4") # high in control and asc


# Tmem69
MEA_dom_DEGs %>% filter(symbol== "Tmem69")# high in control
mPFC_dom_DEGS %>% filter(symbol== "Tmem69")# none
MEA_sub_DEGs %>% filter(symbol== "Tmem69") # none
mPFC_sub_DEGs %>% filter(symbol== "Tmem69") # low in asc

# Sdc4
MEA_dom_DEGs %>% filter(symbol== "Sdc4")# low in doms
mPFC_dom_DEGS %>% filter(symbol== "Sdc4")# low in des
MEA_sub_DEGs %>% filter(symbol== "Sdc4") # none
mPFC_sub_DEGs %>% filter(symbol== "Sdc4") # low in subs

# Cap1
MEA_dom_DEGs %>% filter(symbol== "Cap1")# low in des
mPFC_dom_DEGS %>% filter(symbol== "Cap1")# none
MEA_sub_DEGs %>% filter(symbol== "Cap1") # high in asc
mPFC_sub_DEGs %>% filter(symbol== "Cap1") # high in subs and asc

# Rnf13
MEA_dom_DEGs %>% filter(symbol== "Rnf13")# high in dom
MEA_sub_DEGs %>% filter(symbol== "Rnf13") # none

# Efhd2
MEA_dom_DEGs %>% filter(symbol== "Efhd2")# none
mPFC_dom_DEGS %>% filter(symbol== "Efhd2")# high in dom
MEA_sub_DEGs %>% filter(symbol== "Efhd2") # high in sub
mPFC_sub_DEGs %>% filter(symbol== "Efhd2") # none

# Dusp1
MEA_dom_DEGs %>% filter(symbol== "Dusp1")# high in des
mPFC_dom_DEGS %>% filter(symbol== "Dusp1")# none
MEA_sub_DEGs %>% filter(symbol== "Dusp1") # high in sub
mPFC_sub_DEGs %>% filter(symbol== "Dusp1") # high in sub


# Ier2 - IER gene
MEA_dom_DEGs %>% filter(symbol== "Ier2")# high in des
MEA_sub_DEGs %>% filter(symbol== "Ier2") # high in sub


# Fos
MEA_dom_DEGs %>% filter(symbol== "Fos")# high in des
mPFC_dom_DEGS %>% filter(symbol== "Fos")# high in dom and des
MEA_sub_DEGs %>% filter(symbol== "Fos") # none
mPFC_sub_DEGs %>% filter(symbol== "Fos") # high in sub

# Arc
MEA_dom_DEGs %>% filter(symbol== "Arc")# high in dom and des
mPFC_dom_DEGS %>% filter(symbol== "Arc")# high in dom
MEA_sub_DEGs %>% filter(symbol== "Arc") # high in sub
mPFC_sub_DEGs %>% filter(symbol== "Arc") # high in sub

# Egr1
MEA_dom_DEGs %>% filter(symbol== "Egr1")# none
mPFC_dom_DEGS %>% filter(symbol== "Egr1")#  high in dom
MEA_sub_DEGs %>% filter(symbol== "Egr1") # high in sub
mPFC_sub_DEGs %>% filter(symbol== "Egr1") # high in sub

# Egr2
MEA_dom_DEGs %>% filter(symbol== "Egr2")# hig in dom ahnd des 
mPFC_dom_DEGS %>% filter(symbol== "Egr2")#  none
MEA_sub_DEGs %>% filter(symbol== "Egr2") # none
mPFC_sub_DEGs %>% filter(symbol== "Egr2") # high in sub




# Bdnf
MEA_dom_DEGs %>% filter(symbol== "Bdnf")# high in dom
mPFC_dom_DEGS %>% filter(symbol== "Bdnf")# high in des
MEA_sub_DEGs %>% filter(symbol== "Bdnf") # none
mPFC_sub_DEGs %>% filter(symbol== "Bdnf") # high in control and sub


# Btg2
mPFC_dom_DEGS %>% filter(symbol== "Btg2")# high in des
mPFC_sub_DEGs %>% filter(symbol== "Btg2") # none

# Nr4a1
MEA_dom_DEGs %>% filter(symbol== "Nr4a1")# high in des
mPFC_dom_DEGS %>% filter(symbol== "Nr4a1")# none
MEA_sub_DEGs %>% filter(symbol== "Nr4a1") # high sub
mPFC_sub_DEGs %>% filter(symbol== "Nr4a1") # high  sub


# klf9
mPFC_dom_DEGS %>% filter(symbol== "Klf9")#  high in des
mPFC_sub_DEGs %>% filter(symbol== "Klf9") # none

# fkbp5
MEA_dom_DEGs %>% filter(symbol== "Fkbp5")# none
mPFC_dom_DEGS %>% filter(symbol== "Fkbp5")# high in des
MEA_sub_DEGs %>% filter(symbol== "Fkbp5") # high in aes
mPFC_sub_DEGs %>% filter(symbol== "Fkbp5") # high in control 


# homer1
MEA_dom_DEGs %>% filter(symbol == "Homer1")# high in dom and des 
mPFC_dom_DEGS %>% filter(symbol== "Homer1")# none
MEA_sub_DEGs %>% filter(symbol== "Homer1") # none
mPFC_sub_DEGs %>% filter(symbol== "Homer1") # high in sub  


# fosb
MEA_dom_DEGs %>% filter(symbol == "Fosb")# high in dom and des 
mPFC_dom_DEGS %>% filter(symbol== "Fosb")# none
MEA_sub_DEGs %>% filter(symbol== "Fosb") # none
mPFC_sub_DEGs %>% filter(symbol== "Fosb") # high in sub  


dom_mpfc_exp
dom_mea_exp

sub_mpfc_exp
sub_mea_exp

# just MEA 

dm <- dom_mea_exp 
sm <- sub_mea_exp
smp <- sub_mpfc_exp 
dmp <- dom_mpfc_exp 



a1<- dm %>% filter(symbol == "Fkbp5")
a2 <- sm %>% filter(symbol == "Fkbp5")
a3 <- smp %>% filter(symbol == "Fkbp5")
a4 <- dmp %>% filter(symbol == "Fkbp5")

b1<- dm %>% filter(symbol == "Neurod6")
b2 <- sm %>% filter(symbol == "Neurod6")
b3 <- smp %>% filter(symbol == "Neurod6")
b4 <- dmp %>% filter(symbol == "Neurod6")

c1 <- dm %>% filter(symbol == "Nxph4")
c2 <- sm %>% filter(symbol == "Nxph4")
c3 <- smp %>% filter(symbol == "Nxph4")
c4 <- dmp %>% filter(symbol == "Nxph4")

d1 <- dm %>% filter(symbol == "Tmem69")
d2 <- sm %>% filter(symbol == "Tmem69")
d3 <- smp %>% filter(symbol == "Tmem69")
d4 <- dmp %>% filter(symbol == "Tmem69")

e1 <- dm %>% filter(symbol == "Mpeg1")
e2 <- sm %>% filter(symbol == "Mpeg1")
e3 <- smp %>% filter(symbol == "Mpeg1")
e4 <- dmp %>% filter(symbol == "Mpeg1")

f1 <- dm %>% filter(symbol == "Mybpc1")
f2 <- sm %>% filter(symbol == "Mybpc1")
f3 <- smp %>% filter(symbol == "Mybpc1")
f4 <- dmp %>% filter(symbol == "Mybpc1")

g1 <- dm %>% filter(symbol == "Arhgap36")
g2 <- sm %>% filter(symbol == "Arhgap36")
g3 <- smp %>% filter(symbol == "Arhgap36")
g4 <- dmp %>% filter(symbol == "Arhgap36")

h1 <- dm %>% filter(symbol == "Ttr")
h2 <- sm %>% filter(symbol == "Ttr")
h3 <- smp %>% filter(symbol == "Ttr")
h4 <- dmp %>% filter(symbol == "Ttr")

dp <-a1 %>% rbind(b1,c1,d1,e1,f1,g1,h1)
dpp <- a4 %>% rbind(b4, c4,d4,e4,f4,g4,h4 )
sp <-a2 %>% rbind(b2,c2,d2,e2,f2,g2,h2 )
spp <- a3 %>% rbind(b3,c3,d3,e3,f3,g3,h3 )

sub_m_id 
sub_p_id 
dom_m_id 
dom_p_id 

dp2 <- dp  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(dom_m_id)

dp2$group <- factor(dp2$group, levels = c("CDOM", "DOM", "DES"))

dpp2 <- dpp  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(dom_p_id)

dpp2$group <- factor(dpp2$group, levels = c("CDOM", "DOM", "DES"))

sp2 <- sp  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(sub_m_id)

sp2$group <- factor(sp2$group, levels = c("CSUB", "SUB", "ASC"))

spp2 <- spp %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(sub_p_id)

spp2$group <- factor(spp2$group, levels = c("CSUB", "SUB", "ASC"))

source('functions/geom_boxjitter.R')

p1 <- ggplot(dp2, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  facet_wrap(~symbol,scales="free_y", ncol = 8)+
  ylab("MEA Normalized Expression") +
  theme_bw()+
  theme(legend.position = "none")
p1


# dpp2$symbol <- factor(dpp2$symbol, levels = c('Dmkn', 'Oxt', 'Cd14', 'Il1b', "Btg2", "Klf9", "Mas1"))


p4 <- ggplot(dpp2, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  facet_wrap(~symbol,scales="free_y", ncol = 8)+
  ylab("mPFC Normalized Expression") +
  theme_bw()+
  theme(legend.position = "none")
p4

p2 <- ggplot(sp2, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  facet_wrap(~symbol,scales="free_y", ncol = 9)+
  ylab("MEA Normalized Expression") +
  theme_bw()+
  theme(legend.position = "none")
p2


p3 <- ggplot(spp2, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  scale_fill_manual(values = c("#f7cb44", "#de7065ff", "#13306dff"))+
  facet_wrap(~symbol,scales="free_y", ncol = 8)+
  ylab("mPFC Normalized Expression") +
  theme_bw()+
  theme(legend.position = "none")
p3
p4

both <- gridExtra::grid.arrange(p1,p4,p2,p3, nrow = 4)


ggsave("manuscript/brain/manuscript70/results/results_figures/bothregions_bothgroups_interesting_genes.png",both, height = 10, width =15, dpi=300)



a2 <- sm %>% filter(symbol == "Chrm3")
b2 <- sm %>% filter(symbol == "Lhx1")
c2 <- sm %>% filter(symbol == "Tmem232")
d2 <- sm %>% filter(symbol == "Ogn")
e2 <- sm %>% filter(symbol == "Serpina1c")
f2 <- sm %>% filter(symbol == "Sim1")
g2 <- sm %>% filter(symbol == "Tmem40")
h2 <- sm %>% filter(symbol == "Ggt1")
j2 <- sm %>% filter(symbol == "Pnmt")
k2 <- sm %>% filter(symbol == "Cdkl1")