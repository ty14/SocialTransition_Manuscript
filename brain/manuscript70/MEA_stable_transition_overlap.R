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

#stable data 
limma_list<- readRDS("manuscript/brain/manuscript/results/Won_MeA_data/limma_MEA.RDS")%>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 


as <- limma_list$alphasub

ab<- limma_list$alphasubdom
  
bs <- limma_list$subdomsub


#Transition data 
limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_MEA_ControlDD.RDS") %>% 
  map(~distinct(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol,ensgene, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

cdom <- limma_list$controldom

cdes <- limma_list$controldes 

dd <- limma_list$domdes 


##overlap 
ab_up <- ab %>% filter(logFC > 0.2) %>% arrange(-logFC)

ab_down <- ab %>% filter(logFC < 0.2) %>% arrange(logFC)

as_up <- as %>% filter(logFC > 0.2) %>% arrange(-logFC)

as_down <- as %>% filter(logFC < 0.2) %>% arrange(logFC)

bs_up <- bs %>% filter(logFC > 0.2) %>% arrange(-logFC)

bs_down <- bs %>% filter(logFC < 0.2) %>% arrange(logFC)



cdom_up <- cdom %>% filter(logFC > 0.2) %>% arrange(-logFC)

cdom_down <- cdom %>% filter(logFC < 0.2) %>% arrange(logFC)


dd_up <- dd %>% filter(logFC > 0.2) %>% arrange(-logFC)

dd_down <- dd%>% filter(logFC < 0.2) %>% arrange(logFC)

#one
cdom_up[ab_up$symbol %in% cdom_up$symbol, ]

cdom_down[ab_down$symbol %in% cdom_down$symbol, ]

cdom[ab$symbol %in% cdom$symbol,]

#two 
cdom_up[as_up$symbol %in% cdom_up$symbol, ]

cdom_down[as_down$symbol %in% cdom_down$symbol, ]

cdom[as$symbol %in% cdom$symbol,]

#three
dd_up[ab_up$symbol %in% dd_up$symbol, ] # look here 

dd_down[ab_down$symbol %in% dd_down$symbol, ] # look here 

dd[ab$symbol %in% dd$symbol,]

#four 
dd_up[as_up$symbol %in% dd_up$symbol, ] # look here 

dd_down[as_down$symbol %in% dd_down$symbol, ] # look here 

dd[ab$symbol %in% dd$symbol,]


#five
cdom_up[bs_up$symbol %in% cdom_up$symbol, ] # look here 
cdom_down[bs_down$symbol %in% cdom_down$symbol, ] # look here 

cdom[bs$symbol %in% dd$symbol,]

#six 
dd_up[bs_up$symbol %in% dd_up$symbol, ] # look here 

dd_down[bs_down$symbol %in% dd_down$symbol, ] # look here 

dd[bs$symbol %in% dd$symbol,]



### Top up-regulated genes
head(as_up, 10) # nothing maybe more related to asc - sub 
head(dd_up, 10)

head(ab_up, 20) # related
head(dd_up, 20)

head(bs_up, 10) # nothing 
head(dd_up, 10)


### Down up-regulated genes
head(as_down, 10) #something 
head(dd_down, 10)

head(ab_down, 10) #maybe?
head(dd_down, 10)

head(bs_down, 10)#nothing 
head(dd_down, 10)

chol <- c('Chrm2', 'Chrm4','Gnai1','Gnai2','Chrm1','Chrm3','Chrm5','Gna11', 'Gna14','Gnaq','Ache',
          'Chat','Grk2', 'Grk5', 'Rgs2', 'Rgs4', 'Rgs6', 'Slc18a3', 'Slc5a7', 'Nat1', 'Lhx8', 'Slc10a4', 
          'Gbx1','Chrna2', 'Chrna3', 'Chrna6', 'Chrna7', 'Chrnb4', 'Chrnb3', 'Agrn', 'Chrna1', 'Chrna10', 
          'Chrna4', 'Chrna5', 'Chrna9', 'Chrnb1', 'Chrnb2', 'Chrnd', 'Chrne', 'Chrng', 'Dok7', 'Lrp4', 'Musk',
          'Rapsn')

x <- as_up[as_up$symbol %in% dd_up$symbol, ]
 x[x$symbol %in% chol,]

y <- ab_up[ab_up$symbol %in% dd_up$symbol, ]
y[y$symbol %in% chol,]

xx <- as_down[as_down$symbol %in% dd_down$symbol, ]
xx[xx$symbol %in% chol,]
yy <- ab_down[ab_down$symbol %in% dd_down$symbol, ]
yy[yy$symbol %in% chol,]


bs[bs$symbol %in% chol,]



#AcH overlap all UP reg: Lhx8, Slc10a4, Slc5a7, Gbx1
# ex <- readRDS("manuscript/brain/manuscript/results/Won_MeA_data/limma_vdl_MEA")

ex <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_MeA_ControlDD.RDS")


dm <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

dom_m_id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

a1<- dm %>% filter(symbol == "Lhx8")
a2 <- dm %>% filter(symbol == "Slc5a7")
a3 <- dm %>% filter(symbol == "Slc10a4")
a4 <- dm %>% filter(symbol == "Gbx1")
a5<- dm %>% filter(symbol == "Chrm2")
a6 <- dm %>% filter(symbol == "Chrna4")
a7 <- dm %>% filter(symbol == "Chrna6")
a8 <- dm %>% filter(symbol == "Gna11")
a9<- dm %>% filter(symbol == "Rgs2")
a10 <- dm %>% filter(symbol == "Agrn")
a11 <- dm %>% filter(symbol == "Ache")
a12 <- dm %>% filter(symbol == "Chat")
a13 <- dm %>% filter(symbol == "Chrm1")
a14 <- dm %>% filter(symbol == "Chrna2")
dp1 <- a1 %>%  rbind(a2, a3,a4, a5,a6,a7,a8,a9,a10,a11,a12,a13,a14)

dp2 <- dp1  %>% pivot_longer(cols = 2:21, names_to = "ids") %>% full_join(dom_m_id)

# dp2$group <- factor(dp2$group, levels = c("Alpha", "Subdominant", "Subordinate")) 

dp2$group <- factor(dp2$group, levels = c("CDOM", "DOM", "DES")) 

dp2 <- dp2 %>% na.omit(.)
source('functions/geom_boxjitter.R')

p1 <- ggplot(dp2, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(3))+
  scale_fill_manual(values = viridis::viridis(3))+
  facet_wrap(~symbol,scales="free_y", ncol = 6)+
  ylab("MEA Normalized Expression") +
  xlab("")+
  theme_bw()+
  theme(legend.position = "none")
p1


p2 <- ggplot(dp2, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  facet_wrap(~symbol,scales="free_y", ncol = 7)+
  ylab("MEA Normalized Expression") +
  xlab("")+
  theme_bw()+
  theme(legend.position = "none")
p2



limma_list<- readRDS("manuscript/brain/manuscript/results/limma_MeA_25hrDD.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 


y<- limma_list$status 
y[y$symbol %in% chol,]
# symbol     logFC    P.Value 
# 92    Chat 0.9534004 0.01134173 
# 226  Chrm2 0.4315808 0.03608892 