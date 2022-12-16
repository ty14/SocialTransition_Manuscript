# PCA and upsetting plots 

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
library(UpSetR)
library(workflowr)
library(ComplexUpset)

#MEA

#DOM
limma_list<- readRDS("brain/results/RDS/limma_MeA_CDOM.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= 0.2)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) %>% map(~filter(.,symbol != ""))

cdom <- limma_list$cdom
cdes <- limma_list$cdes 
domdes <- limma_list$domdes 

cdom_up <- cdom %>% filter(logFC > 0.2) %>% arrange(-logFC)

cdom_down <- cdom %>% filter(logFC < 0.2) %>% arrange(logFC)

cdes_up <- cdes %>% filter(logFC > 0.2) %>% arrange(-logFC)

cdes_down <- cdes %>% filter(logFC < 0.2) %>% arrange(logFC)

dd_up <- domdes %>% filter(logFC > 0.2) %>% arrange(-logFC)

dd_down <- domdes%>% filter(logFC < 0.2) %>% arrange(logFC)

#SUB
limma_list<- readRDS("brain/results/RDS/limma_MeA_CSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= 0.2)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) %>% map(~filter(.,symbol != ""))

cs<- limma_list$csub
ca <- limma_list$casc
sa <- limma_list$subasc


cs_up <- cs %>% filter(logFC > 0.2) %>% arrange(-logFC)

cs_down <- cs %>% filter(logFC < 0.2) %>% arrange(logFC)

ca_up <- ca %>% filter(logFC > 0.2) %>% arrange(-logFC)

ca_down <- ca %>% filter(logFC < 0.2) %>% arrange(logFC)

sa_up <- sa %>% filter(logFC > 0.2) %>% arrange(-logFC)

sa_down <- sa %>% filter(logFC < 0.2) %>% arrange(logFC)




DEG <-  readRDS("brain/results/RDS/limma_MeA_CDOM.RDS") %>% 
  map(~distinct(.)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) %>% map(~filter(.,symbol != ""))

DEGx <- DEG$cdom


# example of list input (list of named vectors)
listInput <- list(CDOM= cdom$symbol, CDES=cdes$symbol, DOMDES = domdes$symbol, 
                  CSUB = cs$symbol, CASC = ca$symbol, SUBASC = sa$symbol)


upset(fromList(listInput), nsets = 6, order.by = "freq", keep.order = F)

#just doms 
DEGx <- cdom %>% rbind(cdes, domdes)

listInput <- list(CDOM_up= cdom_up$symbol, CDES_up=cdes_up$symbol, DOMDES_up = dd_up$symbol,
                  CDOM_down= cdom_down$symbol, CDES_down=cdes_down$symbol, DOMDES_down = dd_down$symbol, 
                  CS_up= cs_up$symbol, CA_up=ca_up$symbol, SA_up = sa_up$symbol,
                  CS_down= cs_down$symbol, CA_down=ca_down$symbol, SA_down = sa_down$symbol)


upset(fromList(listInput), nsets = 12, order.by = "freq", keep.order = F)

wrap_label <- function(x) str_wrap(str_replace_all(x, "_", " "), width = 20)

#doms up and down: 
list_input <- list("CDOM up-regulated" = cdom_up$symbol,
                   "CDOM down-regulated" = cdom_down$symbol,
                   "CDES up-regulated" = cdes_up$symbol,
                   "CDES down-regulated" = cdes_down$symbol,
                   "DD up-regulated" = dd_up$symbol,
                   "DD down-regulated" = dd_down$symbol)



asc <- read_csv("brain/results/tables/ascender_MEA_genes.csv")
des <- read_csv("brain/results/tables/descender_MEA_genes.csv")
dom_dis <- read_csv("brain/results/tables/distruption_MEA_genes.csv")
sub_dis <- read_csv("brain/results/tables/distruptionSUB_MEA_genes.csv")

asc_up <- asc %>% filter(reg == "Up")
asc_down <- asc %>% filter(reg == "Down")
des_up <- des %>% filter(reg == "Up")
des_down <- des %>% filter(reg == "Down")
dom_up <- dom_dis %>% filter(reg == "Up")
dom_down <-dom_dis %>% filter(reg == "Down")
sub_up <- sub_dis %>% filter(reg == "Up")
sub_down <-sub_dis %>% filter(reg == "Down")


list_input <- list("ASC up-regulated" = asc_up$symbol,
                   "ASC down-regulated" = asc_down$symbol,
                   "DES up-regulated" = des_up$symbol,
                   "DES down-regulated" = des_down$symbol,
                   "DOM-dis up-regulated" = dom_up$symbol,
                   "DOM-dis down-regulated" = dom_down$symbol, 
                   "SUB-dis up-regulated" = sub_up$symbol,
                   "SUB-dis down-regulated" = sub_down$symbol)

data <- fromList(list_input)
data

inx <- upset(data,nsets = 8, order.by = "freq", keep.order = F)

des_up$symbol[des_up$symbol %in% dom_up$symbol]   #"Tmem42"
des_up$symbol[des_up$symbol %in% sub_down$symbol] #"Zfp174"
sub_up$symbol[sub_up$symbol %in% dom_down$symbol] #"Zc3h14"
sub_up$symbol[sub_up$symbol %in% dom_up$symbol]   #"Mt3"
sub_up$symbol[sub_up$symbol %in% asc_up$symbol]   #"Lrba"
des_down$symbol[des_down$symbol %in% asc_down$symbol] # "Kif21a" "Nap1l5" "Rai2"
des_up$symbol[des_up$symbol %in% asc_up$symbol]       #"Tesc"    "Unc13a"  "Lpl"     "Car12"   "Pitpnm2"



#aggressive genes?
sa_down$symbol[sa_down$symbol %in% dd_up$symbol]
#"Tgfa"    "Irf5"    "Gjd4"    "Bmpr1b"  "Tmem123" "Cog2"    "Carmil1" "Agap1"  

sa_down$symbol[sa_down$symbol %in% cdom_up$symbol]
# "Atp6v1c2" "Lamb3"    "Rreb1"    "Col23a1"  "Doc2b"    "Gda"      "Mr1"      "Nek10"    "Rfx3"     "Fermt3"   "Fam184b"  "Map3k14" 
# "Olfr109"  "Jade3"    "Plek"     "Fxyd5"    "Fam126a"  "Zfp40"    "Olfr393"  "Il17rd"   "Chst11"   "Pced1b"   "Sox9"     "Mast3"   
# "Vmn1r181" "Ppp3ca"   "Ttc8"     "Csad"     "P3h3"     "Shisa7"   "Capn15"   "Vmn2r113" "Efnb3"    "C2cd2l"   "Neurl1b"  "Cdh13"   
# "Fam227a"  "Stxbp2"   "Jph4"     "Slc7a11"  "Adpgk"    "Dennd1b"  "Ssr3"     "Paqr9"    "Olfr539"  "Prkce"    "Pcnx3"    "Mtr"     
# "Tspoap1"  "Mbip"     "Ikzf4"    "Unc79"    "Carf"     "Arhgap1"  "Uvssa"    "Heatr3"   "Alkbh1"   "Alkbh1"   "Slc25a37" "Tab2"    
# "Tmem260"  "Ncstn"    "Zfp512b"  "Prdm10"   "Grid1"    "Sox12"    "Nhsl2"    "Map3k2"   "Fbrs"     "Sipa1l3"  "Plod3"  

#social defeat?
sa_up$symbol[sa_up$symbol %in% dd_down$symbol]
# "Kcnk12" "Zfp850" "Matn2"  "Eif3j2" "Fads3"  "Lins1" 


sa_up$symbol[sa_up$symbol %in% cdes_down$symbol]
# "Kcnk12"   "Ccne1"    "Aurkaip1" "Atf4"     "Cyth2"    "Ndrg3"    "Drg1"     "Tusc2"  