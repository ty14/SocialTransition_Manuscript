library('igraph')
library('network')
library('networkD3')
# library('intergraph')
library('WGCNA')
library(tidyverse)

library(iDINGO)


data(gbm)
colnames(gbm)



dat <- readRDS("manuscript/brain/manuscript70/results/WGCNA/datExpr_MEA_CSUB.RDS")
colnames(dat)

expr <- t(dat) %>% data.frame(.) %>% tibble::rownames_to_column(var = "ensgene")

MM <- readRDS("manuscript/brain/manuscript70/results/WGCNA/MEA_CDOM_WGCNA_MM_GS_all.RDS")
colnames(MM)

hubgenes_df <- read_csv("manuscript/brain/manuscript70/results/wgcna/wgcna_table/MEA_CSUB_hubgene_list.csv")
head(hubgenes_df)

MEs <- readRDS("manuscript/brain/manuscript70/results/WGCNA/MEA_CSUB_MEs.RDS")
head(MEs)


green <- MM %>% filter(module == "green")
#541

green_hub <- hubgenes_df %>%  filter(moduleName == 'green')
#31


x <- expr[expr$ensgene %in% green_hub$ensgene,]



x <- x %>% left_join(., grcm38 %>% dplyr::select(ensgene, symbol)) %>% 
  unique(.) %>% 
  filter(symbol != "NA")


rownames(x)<-NULL



x1 <- x %>% column_to_rownames(., var = "symbol" ) %>%  select(-ensgene)



ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleName") %>%
  pivot_longer(cols = 2:24, names_to = "Module") %>% 
  full_join(coldata) %>% select(SampleName, condition1) %>%   
  filter(condition1 != "DES") %>%
  filter(condition1 != "CDOM") %>%
  filter(condition1 != "DOM") %>% unique(.)

head(ME_df$condition1)

x2 <- t(x1) %>% data.frame(.) %>% 
  rownames_to_column(var = "SampleName") %>% 
  full_join(ME_df) %>% select(-SampleName) 


ifelse(x2$condition1 == "CSUB", 1, x2$condition1) -> x2$condition2
ifelse(x2$condition2 == "SUB", 0, x2$condition2) -> x2$condition2
ifelse(x2$condition2 == "ASC", -1, x2$condition2) -> x2$condition2


g_hub <- x2 %>% select(-condition1) %>% select(x=condition2,1::87 ) %>% filter(x != 1)


 
 g_hub <- x2 %>% select(-condition1) %>% select(1:87) 
 g_hub2 <- t(g_hub)
 
 g_full <- x2 %>% select(-condition1) %>% select(1:540) 
 g_full2 <- t(g_full)
 
 
 
 
 label <- green_hub %>% left_join(., grcm38 %>% dplyr::select(ensgene, symbol)) %>% 
    select(symbol)
 
 
 label <- t(label)
 
 rownames(label) <- NULL
 colnames(label) <- NULL
 
 library(igraph)

 
 # Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
 g <- graph.adjacency(
   as.matrix(as.dist(cor(t(g_hub2), method="pearson"))),
   mode="undirected",
   weighted=TRUE,
   diag=FALSE
 )

 # Simplfy the adjacency object
 g <- igraph::simplify(g, remove.multiple=T, remove.loops=F)
 
 # Colour negative correlation edges as blue
 E(g)[which(E(g)$weight<0)]$color <- "blue"
 
 # Colour positive correlation edges as red
 E(g)[which(E(g)$weight>0)]$color <- "black"
 
 # Convert edge weights to absolute values
 E(g)$weight <- abs(E(g)$weight)
 
 # Change arrow size
 # For directed graphs only
 #E(g)$arrow.size <- 1.0
 
 # Remove edges below absolute Pearson correlation 0.8
 g <- delete_edges(g, E(g)[which(E(g)$weight<0.8)])
 
 # Remove any vertices remaining that have no edges
 g <- delete.vertices(g, degree(g)==0)
 
 # Assign names to the graph vertices (optional)
 V(g)$name <- V(g)$name

 
 # Change shape of graph vertices
 V(g)$shape <- "sphere"
 
 # Change colour of graph vertices
 V(g)$color <- "light green"
 
 # Change colour of vertex frames
 V(g)$vertex.frame.color <- "white"
 
 # Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
 # Multiply scaled vales by a factor of 10
 scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
 vSizes <- (scale01(apply(g_hub, 1, mean)) + 1.0) * 10
 
 # Amplify or decrease the width of the edges
 edgeweights <- E(g)$weight * 3.0
 
 # Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
 mst <- mst(g, algorithm="prim")
 
 # Plot the tree object
 plot(
   mst,
   layout=layout.fruchterman.reingold,
   edge.curved=TRUE,
   vertex.size=vSizes,
   vertex.label.dist= 1,
   vertex.label.color="black",
   asp=FALSE,
   vertex.label.cex= 1,
   edge.width=edgeweights,
   edge.arrow.mode=0,
   main="MEA: CSUB Green")

 
 dev.off() 
 
 