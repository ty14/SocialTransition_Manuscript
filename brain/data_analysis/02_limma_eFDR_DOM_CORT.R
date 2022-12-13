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



# Expression values
dlNorm <-  read.csv("brain/AMY_counts.csv", row.names = 1)
#remove zeros
dlNorm <- dlNorm[apply(dlNorm[], 1, function(x) !all(x==0)),]
#trim sample ids
colnames(dlNorm)[c(1:67)] <- substr(colnames(dlNorm)[c(1:67)], 7, 13)

#Group traits
#Getting metadata ready 
coldata <- read_csv("brain/sample70min_table.csv")
head(coldata)
str(coldata)

#fixing things
coldata$region <- gsub("mPF", "PFC", coldata$region)
coldata$groupEX <- coldata$group


# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))

coldata <-coldata %>% 
  mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))%>%
  mutate(N_Aggrec =  (AggRec70min - min(AggRec70min,na.rm=T))/(max(AggRec70min,na.rm=T)-min(AggRec70min,na.rm=T)))

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


#just get samples I want
coldata$SampleID <- substr(coldata$SampleName, 7, 13)

coldata <- coldata %>%  dplyr::select(-SampleName) %>% 
  filter(Postrank != 3) %>% 
  filter(condition1 != 'ascenders') %>% 
  pivot_wider(
    names_from = region, 
    values_from = region)

row.names <- coldata$SampleID
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
  dplyr::select(SampleID, post_Ncort, N_Aggrec) -> var_info  

row.names <- var_info$SampleID

row.names(var_info) <- row.names #Assigning row names from as sample names  
head(var_info )

dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check

##following Won's code
var_info$post_Ncort -> cort.dl
var_info$N_Aggrec -> agg.dl

design.dl <- model.matrix(~ 0 + cort.dl+agg.dl)
colnames(design.dl) -> mycolnames

v.dl = voom(dge.dl, design.dl, plot = F)
vfit.dl = lmFit(v.dl, design.dl)

efit.dl2 = eBayes(vfit.dl)

p.dl.limma2 = efit.dl2[["p.value"]]
head(p.dl.limma2)

saveRDS(v.dl, "brain/results/RDS/limma_vdl_MEA_DOM_CORT.RDS")



# How many random sampling
R = 5000

set.seed(312)

p.dl.rand = vector('list',length = R)

for( g in 1 : R){
  print(paste("Starting on Permutation", g))
  
  # Randomize the traits
  cort.dl.rand = sample(cort.dl)
  agg.dl.rand = sample(agg.dl)
  
  # Model
  design.dl.rand = model.matrix(~0 + cort.dl.rand + agg.dl.rand)
  colnames(design.dl.rand) <- mycolnames
  
  # Calculate p-values based on randomized traits
  v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
  vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
  
  efit.dl.rand2 = eBayes(vfit.dl.rand)
  
  p.dl.rand[[g]] = efit.dl.rand2[["p.value"]]
  head(p.dl.rand[[g]])
}
q.dl = matrix(0, nrow = nrow(p.dl.limma2), ncol = ncol(p.dl.limma2))

for(h in 1 : R){
  print(paste("Calculating Permutation", h))
  
  temp = p.dl.rand[[h]]
  
  for(c in 1 : 2){
    for(r in 1 : nrow(p.dl.limma2)){
      if(temp[r, c] <= p.dl.limma2[r, c]){
        q.dl[r, c] = q.dl[r, c] + 1
      }
    }
  }
}



q.dl = q.dl / R
colnames(q.dl) = mycolnames
q.dl = as.data.frame(q.dl)

efit.dl2[["p.value"]] <- q.dl
row.names(q.dl) <- NULL
sum(duplicated(row.names(efit.dl2$coefficients)))

saveRDS(q.dl,("brain/results/RDS/limma_vdl_cutoff5_2000_tworand_MEA_DOM_CORT.RDS"))

limma_list <- list()

topTable(efit.dl2, sort.by = "P", n = Inf, coef = 1) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$cort


topTable(efit.dl2, sort.by = "P", n = Inf, coef = 2) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$agg

saveRDS(limma_list,"brain/results/RDS/limma_MEA_DOM_CORT.RDS")

 ##########################################################################################       
        
       # Expression values
        dlNorm <-  read.csv("brain/PFC_counts.csv", row.names = 1)
        #remove zeros
        dlNorm <- dlNorm[apply(dlNorm[], 1, function(x) !all(x==0)),]
        #trim sample ids
        colnames(dlNorm)[c(1:67)] <- substr(colnames(dlNorm)[c(1:67)], 7, 13)
        
        #Group traits
        #Getting metadata ready 
        coldata <- read_csv("brain/sample70min_table.csv")
        head(coldata)
        str(coldata)
        
        coldata <- coldata %>% 
          filter(SampleName != "B1.PFCB12.3.4.trim.sam.counts") 
        
        #fixing things
        coldata$region <- gsub("mPF", "PFC", coldata$region)
        coldata$groupEX <- coldata$group
        
        
        # Normalizing cort data
        # df <- transform(df, N = (N - min(N)) / (max(N) - min(N))
        
        coldata <-coldata %>% 
          mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))%>%
          mutate(N_Aggrec =  (AggRec70min - min(AggRec70min,na.rm=T))/(max(AggRec70min,na.rm=T)-min(AggRec70min,na.rm=T)))
        
        
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
        
        
        #just get samples I want
        coldata$SampleID <- substr(coldata$SampleName, 7, 13)
        
        coldata <- coldata %>%  dplyr::select(-SampleName) %>% 
          filter(Postrank != 3) %>% 
          filter(condition1 != 'ascenders') %>% 
          pivot_wider(
            names_from = region, 
            values_from = region)
        
        row.names <- coldata$SampleID
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
        # 15112   40
        
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
          dplyr::select(SampleID, post_Ncort, N_Aggrec) -> var_info  
        
        row.names <- var_info$SampleID
        
        row.names(var_info) <- row.names #Assigning row names from as sample names  
        head(var_info )
        
        dlNorm<- dlNorm[, rownames(var_info)]
        all(rownames(var_info) == colnames(dlNorm)) #check
        
        ##following Won's code
        var_info$post_Ncort -> cort.dl
        var_info$N_Aggrec -> agg.dl
        
        design.dl <- model.matrix(~ 0 + cort.dl+agg.dl)
        colnames(design.dl) -> mycolnames
        
        v.dl = voom(dge.dl, design.dl, plot = F)
        vfit.dl = lmFit(v.dl, design.dl)
        
        efit.dl2 = eBayes(vfit.dl)
        
        p.dl.limma2 = efit.dl2[["p.value"]]
        head(p.dl.limma2)
        
        saveRDS(v.dl, "brain/results/RDS/limma_vdl_mPFC_DOM_CORT.RDS")
        
        
        
        # How many random sampling
        R = 5000
        
        set.seed(312)
        
        p.dl.rand = vector('list',length = R)
        
        for( g in 1 : R){
          print(paste("Starting on Permutation", g))
          
          # Randomize the traits
          cort.dl.rand = sample(cort.dl)
          agg.dl.rand = sample(agg.dl)
          
          # Model
          design.dl.rand = model.matrix(~0 + cort.dl.rand + agg.dl.rand)
          colnames(design.dl.rand) <- mycolnames
          
          # Calculate p-values based on randomized traits
          v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
          vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
          
          efit.dl.rand2 = eBayes(vfit.dl.rand)
          
          p.dl.rand[[g]] = efit.dl.rand2[["p.value"]]
          head(p.dl.rand[[g]])
        }
        q.dl = matrix(0, nrow = nrow(p.dl.limma2), ncol = ncol(p.dl.limma2))
        
        for(h in 1 : R){
          print(paste("Calculating Permutation", h))
          
          temp = p.dl.rand[[h]]
          
          for(c in 1 : 2){
            for(r in 1 : nrow(p.dl.limma2)){
              if(temp[r, c] <= p.dl.limma2[r, c]){
                q.dl[r, c] = q.dl[r, c] + 1
              }
            }
          }
        }
        
        
        
        q.dl = q.dl / R
        colnames(q.dl) = mycolnames
        q.dl = as.data.frame(q.dl)
        
        efit.dl2[["p.value"]] <- q.dl
        row.names(q.dl) <- NULL
        sum(duplicated(row.names(efit.dl2$coefficients)))
        
        saveRDS(q.dl,("brain/results/RDS/limma_vdl_cutoff5_2000_tworand_mPFC_DOM_CORT.RDS"))
        
        limma_list <- list()
        
        topTable(efit.dl2, sort.by = "P", n = Inf, coef = 1) %>% 
          rownames_to_column('ensgene') %>% 
          left_join(grcm38) %>%
          filter(!is.na(symbol)) %>% 
          dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$cort
        
        
        topTable(efit.dl2, sort.by = "P", n = Inf, coef = 2) %>% 
          rownames_to_column('ensgene') %>% 
          left_join(grcm38) %>%
          filter(!is.na(symbol)) %>% 
          dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$agg
        
        saveRDS(limma_list,"brain/results/RDS/limma_mPFC_DOM_CORT.RDS")
        
        
        #quick look at number of genes
        limma_list %>% map(~filter(., P.Value<0.05)) %>% 
          map(~summarise(.,Up = sum(logFC>0.2),
                         Down = sum(logFC<0.2))) %>% 
          map(~mutate(.,Total = Up + Down))
        
        limma_list %>% map(~hist(.$logFC))
        
       