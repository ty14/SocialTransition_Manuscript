library(ViSEAGO)
library(tidyverse)

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_MEA_ControlDD.RDS") %>% 
  map(~distinct(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

cdom <- limma_list$controldom %>% dplyr::select(entrez)

cdes <- limma_list$controldes %>% dplyr::select(entrez)

domdes <- limma_list$domdes %>% dplyr::select(entrez)

#GET background genes:
x <- cdom %>% rbind(cdes, domdes) 
x$entrez <- as.character(x$entrez) 

background <- x[,1]

# load Differentialy Expressed (DE) gene identifiants from lists

# asc <- read_csv('manuscript/brain/manuscript70/results/tables/ascenders_tranisiton_MEA_genes.csv')
# head(asc)
# 
# #cluster 1
# asc <- asc %>% left_join(., grcm38 %>% dplyr::select(symbol, entrez)) %>% dplyr::select(entrez) 
# asc$entrez <- as.character(asc$entrez)
# asc <- asc[["entrez"]]
# 
# #cluster2 
# asc_down <- asc %>% left_join(., grcm38 %>% dplyr::select(symbol, entrez)) %>% 
#   filter(reg == "Down") %>% dplyr::select(entrez)
# asc_down$entrez <- as.character(asc_down$entrez)
# asc_down <- asc_down[["entrez"]]


des <- read_csv('manuscript/brain/manuscript70/results/tables/descenders_tranisiton_MEA_genes.csv')
head(des)

#cluster 3
des_up <- des %>% left_join(., grcm38 %>% dplyr::select(symbol, entrez)) %>% 
  filter(reg == "Up") %>% dplyr::select(entrez) 
des_up$entrez <- as.character(des_up$entrez)
des_up <- des_up[["entrez"]]

#cluster 4
des_down <- des %>% left_join(., grcm38 %>% dplyr::select(symbol, entrez)) %>% 
  filter(reg == "Down") %>% dplyr::select(entrez)
des_down$entrez <- as.character(des_down$entrez)
des_down <- des_down[["entrez"]]

# connect to Bioconductor
Bioconductor<-ViSEAGO::Bioconductor2GO()

# load GO annotations from Bioconductor
myGENE2GO<-ViSEAGO::annotate(
  "org.Mm.eg.db",
  Bioconductor
)

#summary
myGENE2GO



# create topGOdata for BP for each list of DE genes
# asc1<-ViSEAGO::create_topGOdata(
#   geneSel=asc,
#   allGenes=background,
#   gene2GO=myGENE2GO, 
#   ont="BP",
#   nodeSize=5
# )
# 
# asc2<-ViSEAGO::create_topGOdata(
#   geneSel=asc_down,
#   allGenes=background,
#   gene2GO=myGENE2GO,
#   ont="BP",
#   nodeSize=5
# )

des1<-ViSEAGO::create_topGOdata(
  geneSel=des_up,
  allGenes=background,
  gene2GO=myGENE2GO,
  ont="BP",
  nodeSize=5
)

des2<-ViSEAGO::create_topGOdata(
  geneSel=des_down,
  allGenes=background,
  gene2GO=myGENE2GO,
  ont="BP",
  nodeSize=5
)


# perform topGO tests
# elim_asc1<-topGO::runTest(
#   asc1,
#   algorithm ="elim",
#   statistic = "fisher",
#   cutOff=0.01
# )
# 
# elim_asc2<-topGO::runTest(
#   asc2,
#   algorithm ="elim",
#   statistic = "fisher",
#   cutOff=0.01
# )

elim_des1<-topGO::runTest(
 des1,
  algorithm ="elim",
  statistic = "fisher",
  cutOff=0.01
)

elim_des2<-topGO::runTest(
  des2,
  algorithm ="elim",
  statistic = "fisher",
  cutOff=0.01
)

# merge topGO results
sd_Results<-ViSEAGO::merge_enrich_terms(
  cutoff=0.01,
  Input=list(
    # aes_upx=c(
    #   "asc1",
    #   "elim_asc1"
    # ),
    # aes_downx=c(
    #   "asc2",
    #   "elim_asc2"
    # ),
    des_upx=c(
      "des1",
      "elim_des1"
    ),
    des_downx=c(
      "des2",
      "elim_des2"
    )
  )
)

sd_up_go<-ViSEAGO::merge_enrich_terms(
  cutoff=0.01,
  Input=list(des_upx=c(
    "des1",
    "elim_des1"
  )))


sd_down_go<-ViSEAGO::merge_enrich_terms(
  cutoff=0.01,
  Input=list(des_downx=c(
    "des1",
    "elim_des1"
  )))


sd_up_df <- sd_up_go@data %>% as.data.frame()
head(sd_up_df)
colnames(sd_up_df)

sd_up_df <- cbind(sd_up_df, phenotype= 1)
sd_up <- sd_up_df %>% select(1,2,5,6,9,8)

colnames(sd_up) <- c("Name","Description", "pvalue","enrich", "phenotype", "genes")
head(sd_up)

sd_upfilename <-file.path("manuscript/brain/manuscript70/results/tables/MEA_DOM_up_VISEAGO_cluster_enr_results.txt")
write.table(sd_up,sd_upfilename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)

###descenders 
sd_down_df <- sd_down_go@data %>% as.data.frame()

sd_down_df <- cbind(sd_down_df, phenotype= -1)
sd_down <- sd_down_df %>% select(1,2,5,6,9,8)

colnames(sd_down) <- c("Name","Description", "pvalue","enrich", "phenotype", "genes")
head(sd_down)

sd_downfilename <-file.path("manuscript/brain/manuscript70/results/tables/MEA_DOM_down_VISEAGO_cluster_enr_results.txt")
write.table(sd_down,sd_downfilename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)



# asc_Results<-ViSEAGO::merge_enrich_terms(
#   cutoff=0.01,
#   Input=list(
#     aes_upx=c(
#       "asc1",
#       "elim_asc1"
#     )))

# display a summary
sd_Results
# asc_Results

# show table in interactive mode
ViSEAGO::show_table(sd_Results)
# ViSEAGO::show_table(asc_Results)

# barchart of significant (or not) GO terms by comparison
ViSEAGO::GOcount(sd_Results)
# ViSEAGO::GOcount(asc_Results)

# display intersections
ViSEAGO::Upset(
  sd_Results,
  file="manuscript/brain/manuscript70/results/tables/VISEAGO_descenders_tranisiton.csv"
)

# ViSEAGO::Upset(
  # asc_Results,
  # file="upset.xls"
# )


# create GO_SS-class object
myGOs<-ViSEAGO::build_GO_SS(
  gene2GO=myGENE2GO,
  enrich_GO_terms=sd_Results
)

# asc_GOs<-ViSEAGO::build_GO_SS(
#   gene2GO=myGENE2GO,
#   enrich_GO_terms=asc_Results
# )

# compute Semantic Similarity (SS)
myGOs<-ViSEAGO::compute_SS_distances(
  myGOs,
  distance="Wang"
)

# asc_GOs<-ViSEAGO::compute_SS_distances(
#   asc_GOs,
#   distance="Wang"
# )


#MDSplot
ViSEAGO::MDSplot(myGOs)


# GOterms heatmap with the default parameters
Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
  myGOs,
  showIC=TRUE,
  showGOlabels=TRUE,
  GO.tree=list(
    tree=list(
      distance="Wang",
      aggreg.method="ward.D2"
    ),
    cut=list(
      dynamic=list(
        pamStage=TRUE,
        pamRespectsDendro=TRUE,
        deepSplit=2,
        minClusterSize =2
      )
    )
  ),
  samples.tree=NULL
)


# display the heatmap
ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOterms"
)


# calculate semantic similarites between clusters of GO terms
Wang_clusters_wardD2<-ViSEAGO::compute_SS_distances(
  Wang_clusters_wardD2,
  distance="BMA"
)


# MDSplot
ViSEAGO::MDSplot(
  Wang_clusters_wardD2,
  "GOclusters"
)



#get go terms
sd_go <- sd_Results@data %>% as.data.frame()

colnames(sd_go)

sd_go$GO.ID




