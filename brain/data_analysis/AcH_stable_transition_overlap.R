# libraries 
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
library(glue)

# Looking at MeA data across time groups. 


#70min data 
limma_list<- readRDS("brain/results/RDS/limma_MeA_CDOM.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= 0.2)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) %>% map(~filter(.,symbol != ""))

cdom <- limma_list$cdom
cdes <- limma_list$cdes 
domdes <- limma_list$domdes 


#Won's data 
my_logFC_threshold = 0.15

#stable data 
limma_list<- readRDS("brain/results/Won_MeA_data/limma_MEA.RDS")%>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 


as <- limma_list$alphasub

ab<- limma_list$alphasubdom

bs <- limma_list$subdomsub



#AcH genes
chol <- c('Chrm2', 'Chrm4','Gnai1','Gnai2','Chrm1','Chrm3','Chrm5','Gna11', 'Gna14','Gnaq','Ache',
          'Chat','Grk2', 'Grk5', 'Rgs2', 'Rgs4', 'Rgs6', 'Slc18a3', 'Slc5a7', 'Nat1', 'Lhx8', 'Slc10a4', 
          'Gbx1','Chrna2', 'Chrna3', 'Chrna6', 'Chrna7', 'Chrnb4', 'Chrnb3', 'Agrn', 'Chrna1', 'Chrna10', 
          'Chrna4', 'Chrna5', 'Chrna9', 'Chrnb1', 'Chrnb2', 'Chrnd', 'Chrne', 'Chrng', 'Dok7', 'Lrp4', 'Musk',
          'Rapsn')


# Break up expression 
#stable
ab_up <- ab %>% filter(logFC > 0.2) %>% arrange(-logFC)

ab_down <- ab %>% filter(logFC < 0.2) %>% arrange(logFC)

as_up <- as %>% filter(logFC > 0.2) %>% arrange(-logFC)

as_down <- as %>% filter(logFC < 0.2) %>% arrange(logFC)

bs_up <- bs %>% filter(logFC > 0.2) %>% arrange(-logFC)

bs_down <- bs %>% filter(logFC < 0.2) %>% arrange(logFC)

#70min
cdom_up <- cdom %>% filter(logFC > 0.2) %>% arrange(-logFC)

cdom_down <- cdom %>% filter(logFC < 0.2) %>% arrange(logFC)

cdes_up <- cdes %>% filter(logFC > 0.2) %>% arrange(-logFC)

cdes_down <- cdes %>% filter(logFC < 0.2) %>% arrange(logFC)

dd_up <- domdes %>% filter(logFC > 0.2) %>% arrange(-logFC)

dd_down <- domdes%>% filter(logFC < 0.2) %>% arrange(logFC)



#first stable vs 70min
cdom_up[ab_up$symbol %in% cdom_up$symbol, ] # 9 overlap 
cdom_down[ab_down$symbol %in% cdom_down$symbol, ]# 21 overlap

cdom_up[as_up$symbol %in% cdom_up$symbol, ] #30 overlap
cdom_down[as_down$symbol %in% cdom_down$symbol, ]#33 overlap

 cdom_up[bs_up$symbol %in% cdom_up$symbol, ] #37 overlap 
 cdom_down[bs_down$symbol %in% cdom_down$symbol, ] # 61 overlap  


#stable 
ab[ab$symbol %in% chol,] # 6 #Lhx8, Slc10a4, Slc5a7, Gbx1, Rgs2, Ache
as[as$symbol %in% chol,] # 5 #Lhx8, Slc10a4, Slc5a7, Gnai1, Chrna2
bs[bs$symbol %in% chol,] # 0

#70min
cdom[cdom$symbol %in% chol,] 
# 9 Slc18a3(d), Chrna6(d), Chrm2(d), Gna11, Agrn, Slc5a7(d), Lhx8(d), Chrna4(d), Slc10a4(d)

cdes[cdes$symbol %in% chol,] # 2 Chrna4(des), Rgs2(des)

domdes[domdes$symbol %in% chol,] 
# 10 Chrm2, Slc18a3, Slc10a4, Lhx8, Slc5a7, Gna11(des), Chrm1 (des), Ache, Gbx1, Chat

# Lhx8, Slc10a4, Slc5a7, Gbx1, Rgs2, Ache overlap with stable and transition

# all AcH genes in tranisition data set 
# Chrm2, Slc18a3, Slc10a4, Lhx8, Slc5a7, Gna11(des), Chrm1 (des), Ache, Gbx1, Chat, Agrn, Chrna4, Chrna6, Rgs2

# 70min top genes MeA
head(cdom_up, 10)
head(cdom_down, 10)

head(cdes_up, 10)
head(cdes_down, 10)

head(dd_up, 10)
head(dd_down, 10)



ex <- readRDS("brain/results/RDS/limma_vdl_MeA_CDOM.RDS")
head(ex)
x <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

# Lhx8, Slc10a4, Slc5a7, Gbx1, Rgs2, Ache Chrna2 overlap with stable and transition


# ch2 <- x %>% filter(symbol == "Chrna2")
l <- x %>% filter(symbol == "Lhx8")
g <- x %>% filter(symbol == "Gbx1")
s18 <- x %>% filter(symbol == "Slc5a7")
s10 <- x %>% filter(symbol == "Slc10a4")
a <- x %>% filter(symbol == "Ache")
r <- x %>% filter(symbol == "Rgs2")

ht2 <-l %>% rbind(g,s18,s10,a,r)

id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

colnames(ht2)
ht2 <- ht2 %>% pivot_longer(cols = 2:21, names_to = "ids")

p <- ht2 %>% full_join(id)

p$group <- factor(p$group, levels = c("CDOM", "DOM", "DES"))

source('functions/geom_boxjitter.R')

p1 <- ggplot(p, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  facet_wrap(~symbol,scales="free_y", ncol =7)+
  ylab("70 min expression") +
  theme_bw()+
  theme(legend.position = "none")
p1





ex <- readRDS("brain/results/Won_MeA_data/limma_vdl_MeA.RDS")
head(ex)
x <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

# Lhx8, Slc10a4, Slc5a7, Gbx1, Rgs2, Ache Chrna2 overlap with stable and transition


# ch2 <- x %>% filter(symbol == "Chrna2")
l <- x %>% filter(symbol == "Lhx8")
g <- x %>% filter(symbol == "Gbx1")
s18 <- x %>% filter(symbol == "Slc5a7")
s10 <- x %>% filter(symbol == "Slc10a4")
a <- x %>% filter(symbol == "Ache")
r <- x %>% filter(symbol == "Rgs2")

ht2 <- l%>% rbind(g,s18,s10,a,r)

id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

colnames(ht2)
ht2 <- ht2 %>% pivot_longer(cols = 2:34, names_to = "ids")

p <- ht2 %>% full_join(id)

# p$group <- factor(p$group, levels = c("Alpha", "Subdominant", "Subdorinate"))

source('functions/geom_boxjitter.R')
library(viridis)

p1 <- ggplot(p, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(3))+
  scale_fill_manual(values = viridis::viridis(3))+
  facet_wrap(~symbol,scales="free_y", ncol =7)+
  ylab("Normalized Expression") +
  theme_bw()+
  theme(legend.position = "none")
p1

ggsave("brain/results/img/stable_AcH_genes.png", width =18 , height = 3)


# all AcH genes in tranisition data set 
# Chrm2, Slc18a3, Slc10a4, Lhx8, Slc5a7, Gna11(des), Chrm1 (des), Ache, Gbx1, Chat, Agrn, Chrna4, Chrna6, 


ex <- readRDS("brain/results/RDS/limma_vdl_MeA_CDOM.RDS")
head(ex)
x <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 


ch4 <- x %>% filter(symbol == "Chrna4")
ch6 <- x %>% filter(symbol == "Chrna6")
l <- x %>% filter(symbol == "Lhx8")
g <- x %>% filter(symbol == "Gbx1")
s5 <- x %>% filter(symbol == "Slc5a7")
s10 <- x %>% filter(symbol == "Slc10a4")
s18 <- x %>% filter(symbol == "Slc18a3")
a <- x %>% filter(symbol == "Ache")
r <- x %>% filter(symbol == "Rgs2")
a2 <- x %>% filter(symbol == "Agrn")
c <- x %>%  filter(symbol == "Chat")
ch1 <- x %>%  filter(symbol == "Chrm1")
ch2 <- x %>%  filter(symbol == "Chrm2")
g2 <- x %>% filter(symbol == "Gna11")

ht2 <-l %>% rbind(g,s18,s10,a,r, ch4,ch6, a2, c, ch1, ch2, s5, g2 )

id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

colnames(ht2)
ht2 <- ht2 %>% pivot_longer(cols = 2:21, names_to = "ids")

p <- ht2 %>% full_join(id)

p$group <- factor(p$group, levels = c("CDOM", "DOM", "DES"))

source('functions/geom_boxjitter.R')

p1 <- ggplot(p, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF"))+
  facet_wrap(~symbol,scales="free_y", ncol =4)+
  ylab("70 min expression") +
  theme_bw()+
  theme(legend.position = "none")
p1
