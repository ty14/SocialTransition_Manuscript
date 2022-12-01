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
dlNorm <-as.data.frame(a_count)
row.names(dlNorm) <- nrows

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
# coldata$post_Ncort[is.na(coldata$post_Ncort)] = mean(post_Ncort)

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
row.names(dlNorm) <- nrows


# dlNorm <- dlNorm[!is.na(rowSums(dlNorm)),]

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

dge.dl_dom <- dge.dl[, dge.dl$samples$group %in% c("DOM", "ASC","DES", "SUB")]
dge.dl_dom$samples$group <- droplevels(dge.dl_dom$samples$group)
dge.dl_dom$samples$group
dge.dl<- dge.dl_dom
dge.dl$samples$group

coldata %>% 
  filter(condition1 != "CDOM") %>%
  # filter(condition1 != "DES") %>%
  filter(condition1 != "CSUB") %>%
  dplyr::select(SampleName, condition1) -> var_info  

row.names <- var_info$SampleName

row.names(var_info) <- row.names #Assigning row names from as sample names  
head(var_info)

dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check

##following Won's code
var_info$condition1 %>%
  factor(.,levels = c("DOM","ASC","DES", "SUB")) -> group.dl


design.dl <- model.matrix(~ 0 + group.dl)
colnames(design.dl) -> mycolnames

v.dl = voom(dge.dl, design.dl, plot = F)
vfit.dl = lmFit(v.dl, design.dl)

contrast.matrix <- makeContrasts(group.dlDOM-group.dlASC,
                                 group.dlDOM-group.dlSUB,
                                 # group.dlDOM-group.dlDES, 
                                 group.dlASC-group.dlDES,
                                 # group.dlASC-group.dlSUB,
                                 group.dlSUB-group.dlDES,
                                 levels=design.dl)

vfit.dl2 <- contrasts.fit(vfit.dl, contrast.matrix)

efit.dl2 = eBayes(vfit.dl2)

p.dl.limma2 = efit.dl2[["p.value"]]
head(p.dl.limma2)

 saveRDS(v.dl, "manuscript/brain/manuscript70/results/RDS/limma_vdl_MEA_ReorganizedGroups.RDS")



i = 4

# How many random sampling
R = 5000

p.dl.rand = vector('list',length = R)

for(g in 1 : R){
  print(paste("Starting on Permutation", g))

  # Randomize the traits

  group.dl.rand = sample(group.dl)

  # Model
  design.dl.rand = model.matrix(~0 + group.dl.rand)
  colnames(design.dl.rand) <- mycolnames

  # Calculate p-values based on randomized traits
  v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
  vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)

  vfit.dl.rand2 <- contrasts.fit(vfit.dl.rand, contrast.matrix)

  efit.dl.rand2 = eBayes(vfit.dl.rand2)

  p.dl.rand[[g]] = efit.dl.rand2[["p.value"]]
  head(p.dl.rand[[g]])
}

q.dl = matrix(0, nrow = nrow(p.dl.limma2), ncol = ncol(p.dl.limma2))

for(h in 1 : R){
  print(paste("Calculating Permutation", h))

  temp = p.dl.rand[[h]]

  for(c in 1 : 4){
    for(r in 1 : nrow(p.dl.limma2)){
      if(temp[r, c] <= p.dl.limma2[r, c]){
        q.dl[r, c] = q.dl[r, c] + 1
      }
    }
  }
}

q.dl = q.dl / R
colnames(q.dl) <- mycolnames
q.dl = as.data.frame(q.dl)
row.names(q.dl) <- rownames(dge.dl)


saveRDS(q.dl,("manuscript/brain/manuscript70/results/RDS/limma_vdl_cutoff5_2000_tworand_MEA_70_ReorganizedGroups.RDS"))

q.dl <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_cutoff5_2000_tworand_MEA_70_ReorganizedGroups.RDS")

#

png(filename ="manuscript/brain/manuscript70/results/results_figures/eFDR_hist_MEA_R5000_ReorganizedGroups.png",
    width = 18, height = 17, units = "cm", res = 300)
hist(q.dl-p.dl.limma2, main = "CDOM MEA")
invisible(dev.off())


efit.dl2[["p.value"]] <- q.dl
row.names(q.dl) <- NULL
sum(duplicated(row.names(efit.dl2$coefficients)))

tmp1 <- contrasts.fit(efit.dl2, coef = 1) # DOM_ASC
tmp1 <- eBayes(tmp1)
tmp2 <- contrasts.fit(efit.dl2, coef = 2) #DOM-SUB
tmp2 <- eBayes(tmp2)
tmp3 <- contrasts.fit(efit.dl2, coef = 3) #ASC-DES
tmp3 <- eBayes(tmp3)
tmp4 <- contrasts.fit(efit.dl2, coef = 4) #SUB-DES
tmp4 <- eBayes(tmp4)

limma_list <- list()


topTable(tmp1, sort.by = "P", n = Inf) %>%
  rownames_to_column('ensgene') %>%
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) ->limma_list$domasc



topTable(tmp2, sort.by = "P", n = Inf) %>%
  rownames_to_column('ensgene') %>%
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>%
  dplyr::select(symbol,logFC,P.Value,adj.P.Val)  -> limma_list$domsub


topTable(tmp3, sort.by = "P", n = Inf) %>%
  rownames_to_column('ensgene') %>%
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>%
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$ascdes

topTable(tmp4, sort.by = "P", n = Inf) %>%
  rownames_to_column('ensgene') %>%
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>%
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$subdes


saveRDS(limma_list,"manuscript/brain/manuscript70/results/RDS/limma_MEA_ReorganizedGroups.RDS")


#quick look at number of genes
limma_list %>% map(~filter(., P.Value<0.05)) %>% 
  map(~summarise(.,Up = sum(logFC>0.2),
                 Down = sum(logFC<0.2))) %>% 
  map(~mutate(.,Total = Up + Down))

limma_list %>% map(~hist(.$logFC))

#########################################################################################
#REDO EVERYTHING FOR CSUB groups

# Bring in count data for PFC
a_countdata <- read_csv("manuscript/brain/PFC_counts.csv")
a_countdata$X -> nrows
a_count <- a_countdata[,-1]
dlNorm <-as.data.frame(a_count)
row.names(dlNorm) <- nrows

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
# coldata$post_Ncort[is.na(coldata$post_Ncort)] = mean(post_Ncort)

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
  filter(region != "AMY")

#just get samples I want
row.names <- coldata$SampleName

# coldata <- coldata %>% dplyr::select(-SampleName)

row.names(coldata) <- row.names #Assigning row names from as sample names  
head(coldata)

#check before normalizing 
dlNorm<- dlNorm[, rownames(coldata)]
all(rownames(coldata) == colnames(dlNorm))

#normalize and filter with all groups 
row.names(dlNorm) <- nrows


# dlNorm <- dlNorm[!is.na(rowSums(dlNorm)),]

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

dge.dl_dom <- dge.dl[, dge.dl$samples$group %in% c("DOM", "ASC","DES", "SUB")]
dge.dl_dom$samples$group <- droplevels(dge.dl_dom$samples$group)
dge.dl_dom$samples$group
dge.dl<- dge.dl_dom
dge.dl$samples$group

coldata %>% 
  filter(condition1 != "CDOM") %>%
  # filter(condition1 != "DES") %>%
  filter(condition1 != "CSUB") %>%
  dplyr::select(SampleName, condition1) -> var_info  

row.names <- var_info$SampleName

row.names(var_info) <- row.names #Assigning row names from as sample names  
head(var_info)

dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check

##following Won's code
var_info$condition1 %>%
  factor(.,levels = c("DOM","ASC","DES", "SUB")) -> group.dl


design.dl <- model.matrix(~ 0 + group.dl)
colnames(design.dl) -> mycolnames

v.dl = voom(dge.dl, design.dl, plot = F)
vfit.dl = lmFit(v.dl, design.dl)

contrast.matrix <- makeContrasts(group.dlDOM-group.dlASC,
                                 group.dlDOM-group.dlSUB,
                                 # group.dlDOM-group.dlDES, 
                                 group.dlASC-group.dlDES,
                                 # group.dlASC-group.dlSUB,
                                 group.dlSUB-group.dlDES,
                                 levels=design.dl)

vfit.dl2 <- contrasts.fit(vfit.dl, contrast.matrix)

efit.dl2 = eBayes(vfit.dl2)

p.dl.limma2 = efit.dl2[["p.value"]]
head(p.dl.limma2)

saveRDS(v.dl, "manuscript/brain/manuscript70/results/RDS/limma_vdl_mPFC_ReorganizedGroups.RDS")



i = 4

# How many random sampling
R = 5000

p.dl.rand = vector('list',length = R)

for(g in 1 : R){
  print(paste("Starting on Permutation", g))
  
  # Randomize the traits
  
  group.dl.rand = sample(group.dl)
  
  # Model
  design.dl.rand = model.matrix(~0 + group.dl.rand)
  colnames(design.dl.rand) <- mycolnames
  
  # Calculate p-values based on randomized traits
  v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
  vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
  
  vfit.dl.rand2 <- contrasts.fit(vfit.dl.rand, contrast.matrix)
  
  efit.dl.rand2 = eBayes(vfit.dl.rand2)
  
  p.dl.rand[[g]] = efit.dl.rand2[["p.value"]]
  head(p.dl.rand[[g]])
}

q.dl = matrix(0, nrow = nrow(p.dl.limma2), ncol = ncol(p.dl.limma2))

for(h in 1 : R){
  print(paste("Calculating Permutation", h))
  
  temp = p.dl.rand[[h]]
  
  for(c in 1 : 4){
    for(r in 1 : nrow(p.dl.limma2)){
      if(temp[r, c] <= p.dl.limma2[r, c]){
        q.dl[r, c] = q.dl[r, c] + 1
      }
    }
  }
}

q.dl = q.dl / R
colnames(q.dl) <- mycolnames
q.dl = as.data.frame(q.dl)
row.names(q.dl) <- rownames(dge.dl)


saveRDS(q.dl,("manuscript/brain/manuscript70/results/RDS/limma_vdl_cutoff5_2000_tworand_mPFC_70_ReorganizedGroups.RDS"))

q.dl <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_cutoff5_2000_tworand_mPFC_70_ReorganizedGroups.RDS")

#

png(filename ="manuscript/brain/manuscript70/results/results_figures/eFDR_hist_mPFC_R5000_ReorganizedGroups.png",
    width = 18, height = 17, units = "cm", res = 300)
hist(q.dl-p.dl.limma2, main = "DOM mPFC")
invisible(dev.off())


efit.dl2[["p.value"]] <- q.dl
row.names(q.dl) <- NULL
sum(duplicated(row.names(efit.dl2$coefficients)))

tmp1 <- contrasts.fit(efit.dl2, coef = 1) # DOM_ASC
tmp1 <- eBayes(tmp1)
tmp2 <- contrasts.fit(efit.dl2, coef = 2) #DOM-SUB
tmp2 <- eBayes(tmp2)
tmp3 <- contrasts.fit(efit.dl2, coef = 3) #ASC-DES
tmp3 <- eBayes(tmp3)
tmp4 <- contrasts.fit(efit.dl2, coef = 4) #SUB-DES
tmp4 <- eBayes(tmp4)

limma_list <- list()


topTable(tmp1, sort.by = "P", n = Inf) %>%
  rownames_to_column('ensgene') %>%
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) ->limma_list$domasc



topTable(tmp2, sort.by = "P", n = Inf) %>%
  rownames_to_column('ensgene') %>%
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>%
  dplyr::select(symbol,logFC,P.Value,adj.P.Val)  -> limma_list$domsub


topTable(tmp3, sort.by = "P", n = Inf) %>%
  rownames_to_column('ensgene') %>%
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>%
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$ascdes

topTable(tmp4, sort.by = "P", n = Inf) %>%
  rownames_to_column('ensgene') %>%
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>%
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$subdes


saveRDS(limma_list,"manuscript/brain/manuscript70/results/RDS/limma_mPFC_ReorganizedGroups.RDS")


#quick look at number of genes
limma_list %>% map(~filter(., P.Value<0.05)) %>% 
  map(~summarise(.,Up = sum(logFC>0.2),
                 Down = sum(logFC<0.2))) %>% 
  map(~mutate(.,Total = Up + Down))

limma_list %>% map(~hist(.$logFC))
