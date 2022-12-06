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

# Bring in count data for MEA
a_countdata <- read_csv("manuscript/brain/AMY_counts.csv")
a_countdata$X -> nrows
a_count <- a_countdata[,-1]
row.names(a_count) <- nrows
dlNorm <-as.data.frame(a_count)


#Getting metadata ready 
coldata <- read_csv("manuscript/brain/sample70min_table.csv")
head(coldata)
str(coldata)

#group thing 
coldata$groupEX <- coldata$group
# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))

coldata <-coldata %>% mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))
coldata <- coldata %>%  dplyr::select(-group, -period)
coldata$post_Ncort[is.na(coldata$post_Ncort)] = mean(post_Ncort)

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

#just get samples I want
row.names <- coldata$SampleName

# coldata <- coldata %>% dplyr::select(-SampleName)

row.names(coldata) <- row.names #Assigning row names from as sample names  
head(coldata)

#check before normalizing 
dlNorm<- dlNorm[, rownames(coldata)]
all(rownames(coldata) == colnames(dlNorm))

#normalize and filter with all groups 

dlNorm <- dlNorm[!is.na(rowSums(dlNorm)),]

d = apply(dlNorm, 2, as.numeric)
dim(d)

d0= DGEList(d, group = coldata$condition1)
dim(d0)
rownames(d0) <- rownames(dlNorm)
d0 <- calcNormFactors(d0)

cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
#15184    40

# Now take out groups that you want
#DOMs first 
dge.dl$samples$group

dge.dl_dom <- dge.dl[, dge.dl$samples$group %in% c("CDOM", "DOM", "DES")]
dge.dl_dom$samples$group <- droplevels(dge.dl_dom$samples$group)
dge.dl_dom$samples$group
dge.dl<- dge.dl_dom
dge.dl$samples$group

coldata %>% 
  filter(condition1 != "ASC") %>%
  filter(condition1 != "SUB") %>%
  filter(condition1 != "CSUB") %>%
  dplyr::select(SampleName, condition1, post_Ncort) -> var_info  



row.names <- var_info$SampleName

row.names(var_info) <- row.names #Assigning row names from as sample names  
head(var_info )

dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check

#just postCORT
var_info$post_Ncort -> post_Ncort.dl

design.dl <- model.matrix(~ post_Ncort.dl)
colnames(design.dl) -> mycolnames

v.dl = voom(dge.dl, design.dl, plot = F)
vfit.dl = lmFit(v.dl, design.dl)
efit.dl = eBayes(vfit.dl)

p.dl.limma = efit.dl[["p.value"]]

saveRDS(v.dl,"manuscript/brain/manuscript70/results/RDS/limma_vdl_MEA_DOM_CORT.RDS")



# How many random sampling
R = 5000
# R = 2 # for test
i = 3

# post cort or change in cort?


   p.dl.rand = vector('list',length = R)

    for(g in 1 : R){
      print(paste("Starting on Permutation", g))

      # Randomize the traits

      post_Ncort.dl.rand = sample(post_Ncort.dl)

      # Model
      design.dl.rand = model.matrix(~post_Ncort.dl.rand)
      colnames(design.dl.rand) <- mycolnames

      # Calculate p-values based on randomized traits
      v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
      vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)

      efit.dl.rand = eBayes(vfit.dl.rand)

      p.dl.rand[[g]] = efit.dl.rand[["p.value"]]
      head(p.dl.rand[[g]])
    }

    q.dl = matrix(0, nrow = nrow(p.dl.limma), ncol = ncol(p.dl.limma))

    for(h in 1 : R){
      print(paste("Calculating Permutation", h))

      temp = p.dl.rand[[h]]

      for(c in 1 : 2){
        for(r in 1 : nrow(p.dl.limma)){
          if(temp[r, c] <= p.dl.limma[r, c]){
            q.dl[r, c] = q.dl[r, c] + 1
          }
        }
      }
    }

    q.dl = q.dl / R
    colnames(q.dl) <- mycolnames
    q.dl = as.data.frame(q.dl)
    row.names(q.dl) <- rownames(dge.dl)

saveRDS(q.dl,"manuscript/brain/manuscript70/results/RDS/limma_vdl_cutoff5_2000_tworand_MEA_DOM_CORT")

q.dl <- readRDS('manuscript/brain/manuscript70/results/RDS/limma_vdl_cutoff5_2000_tworand_MEA_DOM_CORT')
  
  
  # png(filename = glue("results_figures/eFDR_hist_{my_region}R{R}_bwchangerate.png"),
  #     width = 18, height = 17, units = "cm", res = 600)
  # hist(q.dl-p.dl.limma, main = glue("{my_region}"))
  # invisible(dev.off())
  # 
  # efit.dl[["p.value"]] <- q.dl
  # row.names(q.dl) <- NULL
  
  topTable(efit.dl, sort.by = "P", n = Inf) %>% 
    rownames_to_column('ensgene') %>% 
    left_join(grcm38) %>%
    filter(!is.na(symbol)) %>% 
    dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_post_Ncort
  
 saveRDS(limma_post_Ncort, "manuscript/brain/manuscript70/results/RDS/limma_PostCORT_DOM_MEA.RDS")

 
 limma_post_Ncort$P.Value <- as.numeric( limma_post_Ncort$P.Value )
 ### Explore
 limma_post_Ncort %>% filter(., P.Value<0.05) %>% 
   summarise(.,Up = sum(logFC>0.2),
                  Down = sum(logFC<0.2)) %>% 
  mutate(.,Total = Up + Down)
 
#######################################################
 #Ascenders,subordinate, and controls
 
 # Bring in count data for MEA
 a_countdata <- read_csv("manuscript/brain/AMY_counts.csv")
 a_countdata$X -> nrows
 a_count <- a_countdata[,-1]
 row.names(a_count) <- nrows
 dlNorm <-as.data.frame(a_count)
 
 
 #Getting metadata ready 
 coldata <- read_csv("manuscript/brain/sample70min_table.csv")
 head(coldata)
 str(coldata)
 
 #group thing 
 coldata$groupEX <- coldata$group
 # Normalizing cort data
 # df <- transform(df, N = (N - min(N)) / (max(N) - min(N))
 
 coldata <-coldata %>% mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))
 coldata <- coldata %>%  dplyr::select(-group, -period)
 coldata$post_Ncort[is.na(coldata$post_Ncort)] = mean(post_Ncort)
 
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
 
 #just get samples I want
 row.names <- coldata$SampleName
 
 # coldata <- coldata %>% dplyr::select(-SampleName)
 
 row.names(coldata) <- row.names #Assigning row names from as sample names  
 head(coldata)
 
 #check before normalizing 
 dlNorm<- dlNorm[, rownames(coldata)]
 all(rownames(coldata) == colnames(dlNorm))
 
 #normalize and filter with all groups 
 
 dlNorm <- dlNorm[!is.na(rowSums(dlNorm)),]
 
 d = apply(dlNorm, 2, as.numeric)
 dim(d)
 
 d0= DGEList(d, group = coldata$condition1)
 dim(d0)
 rownames(d0) <- rownames(dlNorm)
 d0 <- calcNormFactors(d0)
 
 cutoff <- 5
 drop <- which(apply(cpm(d0), 1, max) < cutoff)
 dge.dl <- d0[-drop,]
 dim(dge.dl)
 #15184    40
 
 # Now take out groups that you want
 #DOMs first 
 dge.dl$samples$group
 
 dge.dl_dom <- dge.dl[, dge.dl$samples$group %in% c("CSUB", "SUB", "ASC")]
 dge.dl_dom$samples$group <- droplevels(dge.dl_dom$samples$group)
 dge.dl_dom$samples$group
 dge.dl<- dge.dl_dom
 dge.dl$samples$group
 
 coldata %>% 
   filter(condition1 != "DES") %>%
   filter(condition1 != "DOM") %>%
   filter(condition1 != "CDOM") %>%
   dplyr::select(SampleName, condition1, post_Ncort) -> var_info  
 
 
 row.names <- var_info$SampleName
 
 row.names(var_info) <- row.names #Assigning row names from as sample names  
 head(var_info )
 
 dlNorm<- dlNorm[, rownames(var_info)]
 all(rownames(var_info) == colnames(dlNorm)) #check
 
 #just postCORT
 var_info$post_Ncort -> post_Ncort.dl
 
 design.dl <- model.matrix(~ post_Ncort.dl)
 colnames(design.dl) -> mycolnames
 
 v.dl = voom(dge.dl, design.dl, plot = F)
 vfit.dl = lmFit(v.dl, design.dl)
 efit.dl = eBayes(vfit.dl)
 
 p.dl.limma = efit.dl[["p.value"]]
 
 saveRDS(v.dl,"manuscript/brain/manuscript70/results/RDS/limma_vdl_MEA_SUB_CORT.RDS")
 
 
 
 # How many random sampling
 R = 5000
 # R = 2 # for test
 i = 3
 
 # post cort or change in cort?
 
 
 p.dl.rand = vector('list',length = R)
 
 for(g in 1 : R){
   print(paste("Starting on Permutation", g))
   
   # Randomize the traits
   
   post_Ncort.dl.rand = sample(post_Ncort.dl)
   
   # Model
   design.dl.rand = model.matrix(~post_Ncort.dl.rand)
   colnames(design.dl.rand) <- mycolnames
   
   # Calculate p-values based on randomized traits
   v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
   vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
   
   efit.dl.rand = eBayes(vfit.dl.rand)
   
   p.dl.rand[[g]] = efit.dl.rand[["p.value"]]
   head(p.dl.rand[[g]])
 }
 
 q.dl = matrix(0, nrow = nrow(p.dl.limma), ncol = ncol(p.dl.limma))
 
 for(h in 1 : R){
   print(paste("Calculating Permutation", h))
   
   temp = p.dl.rand[[h]]
   
   for(c in 1 : 2){
     for(r in 1 : nrow(p.dl.limma)){
       if(temp[r, c] <= p.dl.limma[r, c]){
         q.dl[r, c] = q.dl[r, c] + 1
       }
     }
   }
 }
 
 q.dl = q.dl / R
 colnames(q.dl) <- mycolnames
 q.dl = as.data.frame(q.dl)
 row.names(q.dl) <- rownames(dge.dl)
 
 saveRDS(q.dl,"manuscript/brain/manuscript70/results/RDS/limma_vdl_cutoff5_2000_tworand_MEA_SUB_CORT")
 
 q.dl <- readRDS('manuscript/brain/manuscript70/results/RDS/limma_vdl_cutoff5_2000_tworand_MEA_SUB_CORT')
 
 
 # png(filename = glue("results_figures/eFDR_hist_{my_region}R{R}_bwchangerate.png"),
 #     width = 18, height = 17, units = "cm", res = 600)
 # hist(q.dl-p.dl.limma, main = glue("{my_region}"))
 # invisible(dev.off())
 # 
 # efit.dl[["p.value"]] <- q.dl
 # row.names(q.dl) <- NULL
 
 topTable(efit.dl, sort.by = "P", n = Inf) %>% 
   rownames_to_column('ensgene') %>% 
   left_join(grcm38) %>%
   filter(!is.na(symbol)) %>% 
   dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_post_Ncort
 
 saveRDS(limma_post_Ncort, "manuscript/brain/manuscript70/results/RDS/limma_PostCORT_SUB_MEA.RDS")
 
 
 limma_post_Ncort$P.Value <- as.numeric( limma_post_Ncort$P.Value )
 ### Explore
 limma_post_Ncort %>% filter(., P.Value<0.05) %>% 
   summarise(.,Up = sum(logFC>0.2),
             Down = sum(logFC<0.2)) %>% 
   mutate(.,Total = Up + Down)
 
 