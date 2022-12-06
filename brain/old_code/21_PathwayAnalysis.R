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
library(pathview)
library(msigdbr)
library(tidyverse)
grcm38 # mouse genes

my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_mPFC_ControlDD.RDS") %>% 
  map(~distinct(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez, description))) %>% 
  map(~filter(.,!is.na(entrez))) 

cdom <- limma_list$controldom
cdes <- limma_list$controldes 
domdes <- limma_list$domdes 

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_MEA_ControlSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez, description))) %>% 
  map(~filter(.,!is.na(entrez))) 

csub <- limma_list$controlsub
casc <- limma_list$controlasc 
sa <- limma_list$subasc 

library(pathview)
library(gage)
library(gageData)
library(KEGGREST)

kg.mmu=kegg.gsets("mmu")
kg.mmu.eg=kegg.gsets("mmu", id.type="entrez")
head(kg.mmu$kg.sets)
head(kg.mmu.eg$kg.sets,2)

#mmu00140 Steroid hormone biosynthesis none in dominant MEA
#mmu04062 Chemokine signaling pathway
#mmu04014 Ras signaling pathway
#mmu04068 FoxO signaling pathway
#mmu04724 Glutamatergic synapse
#mmu04725 Cholinergic synapse
#mmu04726 Serotonergic synapse
#mmu04727 GABAergic synapse
#mmu04728 Dopaminergic synapse
#mmu04918 Thyroid hormone synthesis
#mmu04935 Growth hormone synthesis, secretion and action


ACH <- kg.mmu.eg$kg.sets$mmu04062 %>% as.data.frame()
colnames(ACH)[1] <- "entrez"
ACH$entrez <- as.integer(ACH$entrez)
selected_genes <- ACH %>% left_join(., grcm38 %>% dplyr::select(symbol, entrez)) %>% select(symbol)



# cat <-msigdbr(species = "Mus musculus", category = "H") 
# name <- unique(cat$gs_name) 
# name
# 
# # HALLMARK_OXIDATIVE_PHOSPHORYLATION
# # HALLMARK_INFLAMMATORY_RESPONSE
# # HALLMARK_TNFA_SIGNALING_VIA_NFKB
# # HALLMARK_APOPTOSIS
# # HALLMARK_ANDROGEN_RESPONSE
# 
# 
msigdbr(species = "Mus musculus", category ="H" ) %>%
  filter(grepl("HALLMARK_OXIDATIVE_PHOSPHORYLATION",gs_name,ignore.case = T)) -> act_sets

act_sets %>%
  .$gs_name %>% unique() -> my_gs_names
my_gs_names

act_sets %>%
  filter(my_gs_names == gs_name) %>%
  .$gene_symbol -> selected_genes


cdom%>% 
  filter(symbol %in% selected_genes) %>% 
  dplyr::select(symbol, logFC, P.Value) %>% 
  mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                      ifelse(logFC>my_logFC_threshold,"CDOM genes","DOM genes"))) %>%
  mutate(Sig = factor(Sig, levels = c("CDOM genes","DOM genes", "N.S."))) %>% 
  mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>% 
  unique() -> df


# cdes%>% 
#   filter(symbol %in% selected_genes$symbol)%>% 
#   dplyr::select(symbol, logFC, P.Value) %>% 
#   mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
#                       ifelse(logFC>my_logFC_threshold,"CDOM genes","DES genes"))) %>%
#   mutate(Sig = factor(Sig, levels = c("CDOM genes","DES genes", "N.S."))) %>% 
#   mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>% 
#   unique() -> df


domdes%>% 
  filter(symbol %in% selected_genes) %>% 
  dplyr::select(symbol, logFC, P.Value) %>% 
  mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                      ifelse(logFC>my_logFC_threshold,"DOM genes","DES genes"))) %>%
  mutate(Sig = factor(Sig, levels = c("DOM genes","DES genes", "N.S."))) %>% 
  mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>% 
  unique() -> df2





keyvals <- ifelse(
  df$logFC < -my_logFC_threshold & df$P.Value<0.05, 'orange',
  ifelse(df$logFC > my_logFC_threshold & df$P.Value<0.05, 'purple4',
         'black'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'purple4'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'orange'] <- 'low'
table(keyvals)


table(keyvals) %>% 
  as.data.frame() %>% 
  mutate(prop = glue::glue("{round(Freq/sum(Freq),3)*100}")) %>% 
  mutate(text = glue::glue("{Freq} genes ({prop}%)")) -> my_texts


if(is.na(my_texts[3,4])){my_texts[3,4] <- "0 genes"}
if(is.na(my_texts[2,4])){my_texts[2,4] <- "0 genes"}

library(EnhancedVolcano)
library(glue)

# png(filename = glue("results_figures/geneset_volcano/EnVol_{my_tissue}_{my_hallmark}.png"),
#     width = 12, height = 11, units = "cm", res = 600)
df %>% 
  filter(P.Value <0.05 & abs(logFC) >0.2) %>% 
  .$symbol -> for_label
EnhancedVolcano(df,
                lab = df$symbol,
                selectLab = for_label,
                x = 'logFC',
                y = 'P.Value',
                title = (""),
                pCutoff = 0.05,
                FCcutoff = 0.0,
                cutoffLineType = 'blank',
                hline = c(0.05),
                hlineCol = c('grey'),
                hlineType = c( 'dotted'),
                hlineWidth = c(1),
                colCustom = keyvals,
                shape = 21,
                pointSize = 2.5,
                labCol = "black",
                labSize = 3.5)+
  annotate("text", x = .5, y = 6,color = "purple4" , size = 3,
           label = glue("Relatively Higher in CDOM \n{my_texts[3,4]}"))+
  annotate("text", x = -.5, y = 6,color = "orange" , size = 3,
           label = glue("Relatively Higher in DOM\n{my_texts[2,4]}"))+
  scale_x_continuous(limits = c(-.75,.75))+
   scale_y_continuous(limits = c(-0.1,8))+
  theme_bw(base_size = 13)+
  labs(color = "",
       y = bquote(~-Log[10]~italic(eFDR)), title = "OXIDATIVE PHOSPHORYLATION")+
  theme(legend.position = "none",
        plot.subtitle = element_blank()) -> temp




##second

keyvals <- ifelse(
  df2$logFC < -my_logFC_threshold & df2$P.Value<0.05, 'orange',
  ifelse(df2$logFC > my_logFC_threshold & df2$P.Value<0.05, 'purple4',
         'black'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'purple4'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'orange'] <- 'low'
table(keyvals)


table(keyvals) %>% 
  as.data.frame() %>% 
  mutate(prop = glue::glue("{round(Freq/sum(Freq),3)*100}")) %>% 
  mutate(text = glue::glue("{Freq} genes ({prop}%)")) -> my_texts


if(is.na(my_texts[3,4])){my_texts[3,4] <- "0 genes"}
if(is.na(my_texts[2,4])){my_texts[2,4] <- "0 genes"}

library(EnhancedVolcano)
library(glue)

# png(filename = glue("results_figures/geneset_volcano/EnVol_{my_tissue}_{my_hallmark}.png"),
#     width = 12, height = 11, units = "cm", res = 600)
df2 %>% 
  filter(P.Value <0.05 & abs(logFC) >0.2) %>% 
  .$symbol -> for_label
EnhancedVolcano(df2,
                lab = df2$symbol,
                selectLab = for_label,
                x = 'logFC',
                y = 'P.Value',
                title = (""),
                pCutoff = 0.05,
                FCcutoff = 0.0,
                cutoffLineType = 'blank',
                hline = c(0.05),
                hlineCol = c('grey'),
                hlineType = c( 'dotted'),
                hlineWidth = c(1),
                colCustom = keyvals,
                shape = 21,
                pointSize = 2.5,
                labCol = "black",
                labSize = 3.5)+
  annotate("text", x = .5, y = 6,color = "purple4" , size = 3,
           label = glue("Relatively Higher in DOM \n{my_texts[3,4]}"))+
  annotate("text", x = -.5, y = 6,color = "orange" , size = 3,
           label = glue("Relatively Higher in DES\n{my_texts[2,4]}"))+
  scale_x_continuous(limits = c(-.75,.75))+
  scale_y_continuous(limits = c(-0.1,8))+
  theme_bw(base_size = 13)+
  labs(color = "",
       y = bquote(~-Log[10]~italic(eFDR)),  title = "OXIDATIVE PHOSPHORYLATION")+
  theme(legend.position = "none",
        plot.subtitle = element_blank()) -> temp1




temp
temp1









##########subs


# HALLMARK_OXIDATIVE_PHOSPHORYLATION
# HALLMARK_INFLAMMATORY_RESPONSE
# HALLMARK_TNFA_SIGNALING_VIA_NFKB
# HALLMARK_APOPTOSIS
# HALLMARK_ANDROGEN_RESPONSE


# msigdbr(species = "Mus musculus", category ="H" ) %>%
#   filter(grepl("HALLMARK_TNFA_SIGNALING_VIA_NFKB",gs_name,ignore.case = T)) -> act_sets
# 
# act_sets %>%
#   .$gs_name %>% unique() -> my_gs_names
# my_gs_names
# 
# act_sets %>%
#   filter(my_gs_names == gs_name) %>%
#   .$gene_symbol -> selected_genes

# 
# #mmu00140 Steroid hormone biosynthesis none in dominant MEA
# #mmu04062 Chemokine signaling pathway
# #mmu04014 Ras signaling pathway
# #mmu04068 FoxO signaling pathway
# #mmu04724 Glutamatergic synapse
# #mmu04725 Cholinergic synapse
# #mmu04726 Serotonergic synapse
# #mmu04727 GABAergic synapse
# #mmu04728 Dopaminergic synapse
# #mmu04918 Thyroid hormone synthesis
# #mmu04935 Growth hormone synthesis, secretion and action
# 
# 
ACH <- kg.mmu.eg$kg.sets$mmu04935 %>% as.data.frame()
colnames(ACH)[1] <- "entrez"
ACH$entrez <- as.integer(ACH$entrez)
selected_genes <- ACH %>% left_join(., grcm38 %>% dplyr::select(symbol, entrez)) %>% select(symbol)



csub%>% 
  filter(symbol %in% selected_genes$symbol) %>% 
  dplyr::select(symbol, logFC, P.Value) %>% 
  mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                      ifelse(logFC>my_logFC_threshold,"CSUB genes","SUB genes"))) %>%
  mutate(Sig = factor(Sig, levels = c("CSUB genes","SUB genes", "N.S."))) %>% 
  mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>% 
  unique() -> df


# casc%>% 
#   filter(symbol %in% selected_genes$symbol)%>% 
#   dplyr::select(symbol, logFC, P.Value) %>% 
#   mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
#                       ifelse(logFC>my_logFC_threshold,"CSUB genes","ASC genes"))) %>%
#   mutate(Sig = factor(Sig, levels = c("CSUB genes","ASC genes", "N.S."))) %>% 
#   mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>% 
#   unique() -> df


sa%>% 
  filter(symbol %in% selected_genes$symbol) %>% 
  dplyr::select(symbol, logFC, P.Value) %>% 
  mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                      ifelse(logFC>my_logFC_threshold,"SUB genes","ASC genes"))) %>%
  mutate(Sig = factor(Sig, levels = c("SUB genes","ASC genes", "N.S."))) %>% 
  mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>% 
  unique() -> df2





keyvals <- ifelse(
  df$logFC < -my_logFC_threshold & df$P.Value<0.05, 'orange',
  ifelse(df$logFC > my_logFC_threshold & df$P.Value<0.05, 'purple4',
         'black'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'purple4'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'orange'] <- 'low'
table(keyvals)


table(keyvals) %>% 
  as.data.frame() %>% 
  mutate(prop = glue::glue("{round(Freq/sum(Freq),3)*100}")) %>% 
  mutate(text = glue::glue("{Freq} genes ({prop}%)")) -> my_texts


if(is.na(my_texts[3,4])){my_texts[3,4] <- "0 genes"}
if(is.na(my_texts[2,4])){my_texts[2,4] <- "0 genes"}

library(EnhancedVolcano)
library(glue)

# png(filename = glue("results_figures/geneset_volcano/EnVol_{my_tissue}_{my_hallmark}.png"),
#     width = 12, height = 11, units = "cm", res = 600)
df %>% 
  filter(P.Value <0.05 & abs(logFC) >0.2) %>% 
  .$symbol -> for_label
EnhancedVolcano(df,
                lab = df$symbol,
                selectLab = for_label,
                x = 'logFC',
                y = 'P.Value',
                title = (""),
                pCutoff = 0.05,
                FCcutoff = 0.0,
                cutoffLineType = 'blank',
                hline = c(0.05),
                hlineCol = c('grey'),
                hlineType = c( 'dotted'),
                hlineWidth = c(1),
                colCustom = keyvals,
                shape = 21,
                pointSize = 2.5,
                labCol = "black",
                labSize = 3.5)+
  annotate("text", x = .5, y = 6,color = "purple4" , size = 3,
           label = glue("Highly expressed in SUB \n{my_texts[2,4]}"))+
  annotate("text", x = -.5, y = 6,color = "orange" , size = 3,
           label = glue("Highly expressed  in ASC \n{my_texts[3,4]}"))+
  scale_x_continuous(limits = c(-.75,.75))+
  scale_y_continuous(limits = c(-0.1,8))+
  theme_bw(base_size = 13)+
  labs(color = "",
       y = bquote(~-Log[10]~italic(eFDR)), title = "Growth Hormone Synthesis and Secretion")+
  theme(legend.position = "none",
        plot.subtitle = element_blank()) -> sub


##second

keyvals <- ifelse(
  df2$logFC < -my_logFC_threshold & df2$P.Value<0.05, 'orange',
  ifelse(df2$logFC > my_logFC_threshold & df2$P.Value<0.05, 'purple4',
         'black'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'purple4'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'orange'] <- 'low'
table(keyvals)


table(keyvals) %>% 
  as.data.frame() %>% 
  mutate(prop = glue::glue("{round(Freq/sum(Freq),3)*100}")) %>% 
  mutate(text = glue::glue("{Freq} genes ({prop}%)")) -> my_texts


if(is.na(my_texts[3,4])){my_texts[3,4] <- "0 genes"}
if(is.na(my_texts[2,4])){my_texts[2,4] <- "0 genes"}

library(EnhancedVolcano)
library(glue)

# png(filename = glue("results_figures/geneset_volcano/EnVol_{my_tissue}_{my_hallmark}.png"),
#     width = 12, height = 11, units = "cm", res = 600)
df2 %>% 
  filter(P.Value <0.05 & abs(logFC) >0.2) %>% 
  .$symbol -> for_label
EnhancedVolcano(df2,
                lab = df2$symbol,
                selectLab = for_label,
                x = 'logFC',
                y = 'P.Value',
                title = (""),
                pCutoff = 0.05,
                FCcutoff = 0.0,
                cutoffLineType = 'blank',
                hline = c(0.05),
                hlineCol = c('grey'),
                hlineType = c( 'dotted'),
                hlineWidth = c(1),
                colCustom = keyvals,
                shape = 21,
                pointSize = 2.5,
                labCol = "black",
                labSize = 3.5)+
  annotate("text", x = .5, y = 6,color = "purple4" , size = 3,
           label = glue("Highly expressed in CSUB \n{my_texts[3,4]}"))+
  annotate("text", x = -.5, y = 6,color = "orange" , size = 3,
           label = glue("Highly expressed  in SUB \n{my_texts[2,4]}"))+
  scale_x_continuous(limits = c(-.75,.75))+
  scale_y_continuous(limits = c(-0.1,8))+
  theme_bw(base_size = 13)+
  labs(color = "",
       y = bquote(~-Log[10]~italic(eFDR)), title = "Growth Hormone Synthesis and Secretion")+
  theme(legend.position = "none",
        plot.subtitle = element_blank()) -> sub2


sub
sub2

