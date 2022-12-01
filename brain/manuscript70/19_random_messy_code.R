library(annotables)
library(clusterProfiler)
library(tidyverse)
grcm38 # mouse genes




datExpr <- readRDS("manuscript/brain/manuscript70/results/wgcna/datExpr_MEA_CDOM.RDS")
net <- readRDS("manuscript/brain/manuscript70/results/wgcna/net_MEA_CDOM_Power9.RDS")
MEs <- readRDS("manuscript/brain/manuscript70/results/wgcna/MEA_CDOM_MEs.RDS")




coldata <- read_csv("manuscript/brain/sample70min_table.csv")
head(coldata)
str(coldata)



coldata$groupEX <- coldata$group
# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))
coldata <-coldata %>% mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))
coldata <- coldata %>%  dplyr::select(-group, -period)
#getting condition1
table(coldata$condition)
# CDOM, RDOM to Descenders (DOM to SUB)(4->1)
# CSUB, SUB to Ascenders (Sub to DOM)  (1->4)
coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "DOM", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "SUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "DES", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "ASC", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 4, "CSUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 1, "CDOM", coldata$condition1)

coldata <- coldata %>%
  filter(Postrank != 3) %>%
  filter(condition1 != 'ascenders') %>%
  filter(region == "AMY")

ME_dom <-MEs%>% data.frame() %>%
  tibble::rownames_to_column(var = "SampleName") %>%
  dplyr::select(SampleName, MEbrown, MEpink, MEpurple, MEgreen) %>%
  pivot_longer(cols = 2:5, names_to = "Module") %>%
  full_join(coldata)
head(ME_dom)

ME_dom$condition1
ME_dom$Module <- gsub("ME", "", ME_dom$Module)

ME_dom  <- ME_dom %>%
  mutate(status = condition1) %>%
  filter(condition1 != "ASC") %>%
  filter(condition1 != "SUB") %>%
  filter(condition1 != "CSUB")

ME_dom$status <- factor(ME_dom$status, levels = c("CDOM", "DOM", "DES"))

modNames = ME_dom$Module

moduleColors = ME_dom$Module

module = "green"
module = "pink"
module = "brown"
module = "purple"
module_list = c( "pink","brown", "green", "purple")
## ==============================================
gettop10GO_WGCNA <- function(module,my_showCategory = 10){
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  colnames(datExpr)[moduleColors==module] -> module_gene
  grcm38 %>%
    filter(ensgene %in% module_gene) %>%
    filter(!is.na(entrez)) %>%
    dplyr::select(entrez) -> go_df_wgcna
  ggo <- enrichGO(gene = go_df_wgcna$entrez %>% unique(),
                  OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                  keyType = "ENTREZID",
                  ont = 'BP',
                  readable = T,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.50)
  fortify(
    ggo,
    showCategory = my_showCategory,
    by = "Count",
    split = NULL,
    includeAll = TRUE
  ) %>%
    dplyr::arrange(desc(GeneRatio)) %>%
    mutate(module = module) -> temp1
  return(rbind(temp1))
}
my_ont = "BP"
my_showCategory = 100
moduleColors %>% unique() -> allcolors
WGCNA_GOs <- vector('list', length(allcolors))
for(i in 1:length(allcolors)){
  gettop10GO_WGCNA(allcolors[i],my_showCategory) -> WGCNA_GOs[[i]]
}


WGCNA_GOs %>%
  do.call(rbind,.) -> MEA_doms_wgcna_all_gos






datExpr <- readRDS("manuscript/brain/manuscript70/results/wgcna/datExpr_MEA_CSUB.RDS")
net <- readRDS("manuscript/brain/manuscript70/results/wgcna/net_MEA_CSUB_Power10.RDS")
MEs <- readRDS("manuscript/brain/manuscript70/results/wgcna/MEA_CSUB_MEs.RDS")




coldata <- read_csv("manuscript/brain/sample70min_table.csv")
head(coldata)
str(coldata)



coldata$groupEX <- coldata$group
# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))
coldata <-coldata %>% mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))
coldata <- coldata %>%  dplyr::select(-group, -period)
#getting condition1
table(coldata$condition)
# CDOM, RDOM to Descenders (DOM to SUB)(4->1)
# CSUB, SUB to Ascenders (Sub to DOM)  (1->4)
coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "DOM", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "SUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "DES", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "ASC", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 4, "CSUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 1, "CDOM", coldata$condition1)

coldata <- coldata %>%
  filter(Postrank != 3) %>%
  filter(condition1 != 'ascenders') %>%
  filter(region == "AMY")

ME_dom <-MEs%>% data.frame() %>%
  tibble::rownames_to_column(var = "SampleName") %>%
  dplyr::select(SampleName, MEyellow, MEsalmon, MEturquoise, MEgreen) %>%
  pivot_longer(cols = 2:5, names_to = "Module") %>%
  full_join(coldata)
head(ME_dom)

ME_dom$condition1
ME_dom$Module <- gsub("ME", "", ME_dom$Module)

ME_dom  <- ME_dom %>%
  mutate(status = condition1) %>%
  filter(condition1 != "DES") %>%
  filter(condition1 != "DOM") %>%
  filter(condition1 != "CDOM")

ME_dom$status <- factor(ME_dom$status, levels = c("CSUB", "DOM", "SUB"))

modNames = ME_dom$Module

moduleColors = ME_dom$Module

module = "yellow"
module = "salmon"
module = "turquoise"
module = "green"
module_list = c( "yellow","salmon", "turquoise", "green")
## ==============================================
gettop10GO_WGCNA <- function(module,my_showCategory = 10){
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  colnames(datExpr)[moduleColors==module] -> module_gene
  grcm38 %>%
    filter(ensgene %in% module_gene) %>%
    filter(!is.na(entrez)) %>%
    dplyr::select(entrez) -> go_df_wgcna
  ggo <- enrichGO(gene = go_df_wgcna$entrez %>% unique(),
                  OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                  keyType = "ENTREZID",
                  ont = 'BP',
                  readable = T,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.50)
  fortify(
    ggo,
    showCategory = my_showCategory,
    by = "Count",
    split = NULL,
    includeAll = TRUE
  ) %>%
    dplyr::arrange(desc(GeneRatio)) %>%
    mutate(module = module) -> temp1
  return(rbind(temp1))
}
my_ont = "BP"
my_showCategory = 100
moduleColors %>% unique() -> allcolors
WGCNA_GOs <- vector('list', length(allcolors))
for(i in 1:length(allcolors)){
  gettop10GO_WGCNA(allcolors[i],my_showCategory) -> WGCNA_GOs[[i]]
}


WGCNA_GOs %>%
  do.call(rbind,.) -> MEA_subs_wgcna_all_gos






#FIRST PINK - salmon
sub <- MEA_subs_wgcna_all_gos %>% dplyr::select(Description, Count, module) %>% 
  filter(module == "yellow") %>% mutate(condition = "SUB") %>% head(5)

# colnames(sub)[2] <- "sub_count"
# colnames(sub)[3] <- "sub_module"


dom <- MEA_doms_wgcna_all_gos %>% dplyr::select(Description, Count, module) %>% 
  filter(module == "brown") %>% mutate(condition = "DOM") %>% head(5)

# colnames(dom)[2] <- "dom_count"
# colnames(dom)[3] <- "dom_module"

all <- sub %>%  rbind(dom)


 p1 <- ggplot(all, aes(x = reorder(Description, -Count), y = Count, fill = condition ) ) +
  geom_col(position = position_dodge(.60), alpha = 2, width = 0.5) +
  labs(y = "Gene Number",
       x = "")+
  scale_color_manual(values = c("dark golden rod", "pale golden rod" ))+
  scale_fill_manual(values = c("dark golden rod", "pale golden rod"))+
  theme_classic()+
  coord_flip()+
  theme(text = element_text(size = 15),
        legend.position = "none")

p1

 ggsave("manuscript/brain/manuscript70/results/results_figures/brown_yellow_MEA_preservation_GO.png",p1, height = 5, width = 8, dpi=300)
 
 
 


 #FIRST PINK - salmon
 sub <- MEA_subs_wgcna_all_gos %>% dplyr::select(Description, Count, module) %>% 
   filter(module == "salmon") %>% mutate(condition = "SUB") %>% head(5)
 
 # colnames(sub)[2] <- "sub_count"
 # colnames(sub)[3] <- "sub_module"
 
 
 dom <- MEA_doms_wgcna_all_gos %>% dplyr::select(Description, Count, module) %>% 
   filter(module == "pink") %>% mutate(condition = "DOM") %>% head(5)
 
 # colnames(dom)[2] <- "dom_count"
 # colnames(dom)[3] <- "dom_module"
 
 all <- sub %>%  rbind(dom)
 
 
 p1 <- ggplot(all, aes(x = reorder(Description, -Count), y = Count, fill = condition ) ) +
   geom_col(position = position_dodge(.60), alpha = 2, width = 0.5) +
   labs(y = "Gene Number",
        x = "")+
   scale_color_manual(values = c("pink", "salmon" ))+
   scale_fill_manual(values = c("pink", "salmon"))+
   theme_classic()+
   coord_flip()+
   theme(text = element_text(size = 15),
         legend.position = "none")
 
 
 
 ggsave("manuscript/brain/manuscript70/results/results_figures/pink_salmon_MEA_preservation_GO.png",p1, height = 5, width = 8, dpi=300)
 
 
 # purple - turquoise
 
 sub <- MEA_subs_wgcna_all_gos %>% dplyr::select(Description, Count, module) %>% 
   filter(module == "turquoise") %>% mutate(condition = "SUB") %>% head(5)
 
 # colnames(sub)[2] <- "sub_count"
 # colnames(sub)[3] <- "sub_module"
 
 
 dom <- MEA_doms_wgcna_all_gos %>% dplyr::select(Description, Count, module) %>% 
   filter(module == "purple") %>% mutate(condition = "DOM") %>% head(5)
 
 # colnames(dom)[2] <- "dom_count"
 # colnames(dom)[3] <- "dom_module"
 
 all <- sub %>%  rbind(dom)
 
 
 p1 <- ggplot(all, aes(x = reorder(Description, -Count), y = Count, fill = condition ) ) +
   geom_col(position = position_dodge(.60), alpha = .5, width = 0.5) +
   labs(y = "Gene Number",
        x = "")+
   scale_color_manual(values = c("purple", "turquoise" ))+
   scale_fill_manual(values = c("purple", "turquoise"))+
   theme_classic()+
   coord_flip()+
   theme(text = element_text(size = 15),
         legend.position = "none")
 
 
 p1
 ggsave("manuscript/brain/manuscript70/results/results_figures/purple_turq_MEA_preservation_GO.png",p1, height = 5, width = 10, dpi=300)
 

 
 # purple - turquoise
 
 sub <- MEA_subs_wgcna_all_gos %>% dplyr::select(Description, Count, module) %>% 
   filter(module == "green") %>% mutate(condition = "SUB") %>% head(5)
 
 # colnames(sub)[2] <- "sub_count"
 # colnames(sub)[3] <- "sub_module"
 
 
 dom <- MEA_doms_wgcna_all_gos %>% dplyr::select(Description, Count, module) %>% 
   filter(module == "green") %>% mutate(condition = "DOM") %>% head(5)
 
 # colnames(dom)[2] <- "dom_count"
 # colnames(dom)[3] <- "dom_module"
 
 all <- sub %>%  rbind(dom)
 
 
 p1 <- ggplot(all, aes(x = reorder(Description, -Count), y = Count, fill = condition ) ) +
   geom_col(position = position_dodge(.60), alpha = .5, width = 0.5) +
   labs(y = "Gene Number",
        x = "")+
   scale_color_manual(values = c("dark sea green", "forestgreen" ))+
   scale_fill_manual(values = c("dark sea green", "forestgreen"))+
   theme_classic()+
   coord_flip()+
   theme(text = element_text(size = 15),
         legend.position = "none")
 
 
 p1
 ggsave("manuscript/brain/manuscript70/results/results_figures/green_green_MEA_preservation_GO.png",p1, height = 5, width = 8, dpi=300)
 
 
 
 limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_MEA_ControlDD.RDS") %>% 
   map(~distinct(.)) %>% 
   map(~filter(.,abs(logFC) >= 0.2)) %>%
   map(~filter(.,P.Value <0.05)) %>% 
   map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
   map(~filter(.,!is.na(entrez))) 
 
 y1a <- limma_list$controldom %>% mutate(contrast = "control-dom")
 
 y2a <- limma_list$controldes %>% mutate(contrast = "control-des")
 
 y3a <- limma_list$domdes %>% mutate(contrast = "dom - des")
 
 MEA_dom_DEGs <- y1a %>% rbind(y2a, y3a)

 MDD <-  MEA_dom_DEGs %>% dplyr::select(symbol) %>% unique(.)
 
 limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_MEA_Controlsub.RDS") %>% 
   map(~distinct(.)) %>% 
   map(~filter(.,abs(logFC) >= 0.2)) %>%
   map(~filter(.,P.Value <0.05)) %>% 
   map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
   map(~filter(.,!is.na(entrez))) 
 
 y1b <- limma_list$controlsub %>% mutate(contrast = "control-sub")
 
 y2b <- limma_list$controlasc %>% mutate(contrast = "control-asc")
 
 y3b <- limma_list$subasc %>% mutate(contrast = "sub - asc")
 
 MEA_sub_DEGs <- y1b %>% rbind(y2b, y3b)
 
 MDS <-  MEA_sub_DEGs %>% dplyr::select(symbol) %>% unique(.)
 
 
brown_syn <-  MEA_doms_wgcna_all_gos %>% filter(Description =="synapse organization") %>% filter(module == "brown")
 syn <- as.character(brown_syn$geneID) %>% strsplit(syn, split = "/")
syn <- do.call(cbind, syn) %>% as.data.frame() %>% select(symbol = V1) 
syn$symbol[(syn$symbol %in% MDD$symbol)] %>% as.data.frame(.)
#brown synapse organization has 42 DEGs for DOMs/105

yellow_syn <-  MEA_subs_wgcna_all_gos %>% filter(Description =="synapse organization") %>% filter(module == "yellow")
syel <- as.character(yellow_syn$geneID) %>% strsplit(syn, split = "/")
syel <- do.call(cbind, syel) %>% as.data.frame() %>% select(symbol = V1) 

syel$symbol[(syel$symbol %in% MDS$symbol)] %>% as.data.frame(.)
#yellow synapse organization has 26 DEGs for DOMs/105


brown_syn <-  MEA_doms_wgcna_all_gos %>% filter(Description =="axonogenesis") %>% filter(module == "brown")
syn <- as.character(brown_syn$geneID) %>% strsplit(syn, split = "/")
syn <- do.call(cbind, syn) %>% as.data.frame() %>% select(symbol = V1) 

syn$symbol[(syn$symbol %in% MDD$symbol)] %>% as.data.frame(.)
#brown anonogensis has 42 DEGs for DOMs/97

yellow_syn <-  MEA_subs_wgcna_all_gos %>% filter(Description =="axonogenesis") %>% filter(module == "yellow")
syel <- as.character(yellow_syn$geneID) %>% strsplit(syn, split = "/")
syel <- do.call(cbind, syel) %>% as.data.frame() %>% select(symbol = V1) 

syel$symbol[(syel$symbol %in% MDS$symbol)] %>% as.data.frame(.)
#yellow synapse organization has 18 DEGs for DOMs/97


brown_syn <-  MEA_doms_wgcna_all_gos %>% filter(Description =="regulation of membrane potential") %>% filter(module == "brown")
syn <- as.character(brown_syn$geneID) %>% strsplit(syn, split = "/")
syn <- do.call(cbind, syn) %>% as.data.frame() %>% select(symbol = V1) 

syn$symbol[(syn$symbol %in% MDD$symbol)] %>% as.data.frame(.)
#brown anonogensis has 35 DEGs for DOMs/ 92


yellow_syn <-  MEA_subs_wgcna_all_gos %>% filter(Description =="regulation of membrane potential") %>% filter(module == "yellow")
syel <- as.character(yellow_syn$geneID) %>% strsplit(syn, split = "/")
syel <- do.call(cbind, syel) %>% as.data.frame() %>% select(symbol = V1) 

syel$symbol[(syel$symbol %in% MDS$symbol)] %>% as.data.frame(.)
#yellow regulation of membrane potential has 21 DEGs for DOMs/97

brown_syn <-  MEA_doms_wgcna_all_gos %>% filter(Description =="regulation of cell growth") %>% filter(module == "brown")
syn <- as.character(brown_syn$geneID) %>% strsplit(syn, split = "/")
syn <- do.call(cbind, syn) %>% as.data.frame() %>% select(symbol = V1) 

syn$symbol[(syn$symbol %in% MDD$symbol)] %>% as.data.frame(.)
#brown anonogensis has  18 DEGs for DOMs/ 91


yellow_syn <-  MEA_subs_wgcna_all_gos %>% filter(Description =="regulation of cell growth") %>% filter(module == "yellow")
syel <- as.character(yellow_syn$geneID) %>% strsplit(syn, split = "/")
syel <- do.call(cbind, syel) %>% as.data.frame() %>% select(symbol = V1) 

syel$symbol[(syel$symbol %in% MDS$symbol)] %>% as.data.frame(.)
#yellow rregulation of cell growth has 9 DEGs for DOMs/97


brown_syn <-  MEA_doms_wgcna_all_gos %>% filter(Description =="positive regulation of cellular catabolic process") %>% filter(module == "brown")
syn <- as.character(brown_syn$geneID) %>% strsplit(syn, split = "/")
syn <- do.call(cbind, syn) %>% as.data.frame() %>% select(symbol = V1) 

syn$symbol[(syn$symbol %in% MDD$symbol)] %>% as.data.frame(.)
#brown positive regulation of cellular catabolic process has 25 DEGs for DOMs/ 90


yellow_syn <-  MEA_subs_wgcna_all_gos %>% filter(Description =="positive regulation of cellular catabolic process") %>% filter(module == "yellow")
syel <- as.character(yellow_syn$geneID) %>% strsplit(syn, split = "/")
syel <- do.call(cbind, syel) %>% as.data.frame() %>% select(symbol = V1) 

syel$symbol[(syel$symbol %in% MDS$symbol)] %>% as.data.frame(.)
#yellow positive regulation of cellular catabolic process has 16 DEGs for DOMs/97

##############################################
  
  sub <- MEA_subs_wgcna_all_gos  %>% 
  filter(module == "salmon") %>% mutate(condition = "SUB") 

# colnames(sub)[2] <- "sub_count"
# colnames(sub)[3] <- "sub_module"


dom <- MEA_doms_wgcna_all_gos %>% 
  filter(module == "pink") %>% mutate(condition = "DOM") 

# colnames(dom)[2] <- "dom_count"
# colnames(dom)[3] <- "dom_module"

all <- sub %>%  rbind(dom) %>% arrange(-Count)

all2 <- head(all,20)


Dpfc <- sub %>% select(Description,GeneRatio,p.adjust,condition) %>% arrange(p.adjust) 
Dpfc <- head(Dpfc)
Damy <- dom %>% select(Description,GeneRatio,p.adjust,condition) %>% arrange(p.adjust)
Damy <- head(Damy)
DD <- Damy %>% full_join(Dpfc)


S1<- ggplot(DD, aes(x= condition, y = fct_reorder(Description, GeneRatio), size=GeneRatio, color=p.adjust, group=condition)) + geom_point(alpha = 0.8) + 
  theme_classic()+
  xlab("") + ylab("")+ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))+
  viridis::scale_color_viridis(option="viridis")

S1 



kin <- readRDS("manuscript/brain/manuscript70/results/kIN/MEA_CDOM_kIN_dataframe_.RDS")

kin_df <- kin %>%  
  arrange(desc(kWithin)) %>% 
  mutate(kINrank = rank(kWithin)) %>% filter(module == "pink")


round(nrow(kin_df)*0.20) -> cutoff

kin_df %>% 
  head(cutoff) %>% 
  .$symbol -> brown_dom



kin <- readRDS("manuscript/brain/manuscript70/results/kIN/MEA_CSUB_kIN_dataframe_.RDS")

kin_df <- kin %>%  
  arrange(desc(kWithin)) %>% 
  mutate(kINrank = rank(kWithin)) %>% filter(module == "salmon")


round(nrow(kin_df)*0.20) -> cutoff

kin_df %>% 
  head(cutoff) %>% 
  .$symbol -> yellow_sub



xx <- brown_dom[brown_dom%in%yellow_sub]


MDD$symbol[MDD$symbol %in% xx]

