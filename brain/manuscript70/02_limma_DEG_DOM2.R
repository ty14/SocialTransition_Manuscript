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


#MEA
a_countdata <- read_csv("manuscript/brain/AMY_counts.csv")
colnames(a_countdata)[1] <- "ensgene"
colnames(a_countdata)[c(2:68)] <- substr(colnames(a_countdata)[c(2:68)], 7, 13)
# #mPFC
# p_countdata <- read_csv("manuscript/brain/PFC_counts.csv")
# colnames(p_countdata)[1] <- "ensgene"
# colnames(p_countdata)[c(2:68)] <- substr(colnames(p_countdata)[c(2:68)], 7, 13)

dlNorm_list <- list(a_countdata )
names(dlNorm_list) <- c("MeA")

#Getting metadata ready 
coldata <- read_csv("manuscript/brain/sample70min_table.csv")
head(coldata)
str(coldata)

#fixing things
coldata$region <- gsub("mPF", "PFC", coldata$region)
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
  filter(SampleName != "B1.PFCB12.3.4.trim.sam.counts") 

#just get samples I want
coldata$SampleID <- substr(coldata$SampleName, 7, 13)

coldata <- coldata %>%  select(-SampleName) %>% 
  pivot_wider(
    names_from = region, 
    values_from = region)

row.names <- coldata$SampleID
row.names(coldata) <- row.names #Assigning row names from as sample names  
head(coldata)


# dlNorm_a <- a_countdata 
# a_countdata<- a_countdata[, rownames(coldata)]
# all(rownames(coldata) == colnames(a_countdata)) #check
# 
# dlNorm_p <- p_countdata 
## p_countdata<- p_countdata[, rownames(coldata)]
# all(rownames(coldata) == colnames(p_countdata)) #check 

# dlNorm_list <- list(dlNorm_a , dlNorm_p)
# names(dlNorm_list) <- c("MeA", "mPFC")

regions <- c("MeA")  

i = 1

# How many random sampling
R = 5000


for (i in 1:length(regions)){
  
  regions[[i]] -> my_region
  
  
  dlNorm <- dlNorm_list[[i]] %>% 
    column_to_rownames('ensgene')
  
  
  colnames(dlNorm)
  
  
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
  
  
  dge.dl$samples$group
  
  dge.dl_dom <- dge.dl[, dge.dl$samples$group %in% c("CDOM", "DOM", "DES")]
  dge.dl_dom$samples$group <- droplevels(dge.dl_dom$samples$group)
  dge.dl_dom$samples$group
  dge.dl<- dge.dl_dom
  dge.dl$samples$group
  
  var_info <- coldata %>% 
    filter(condition1 != "SUB") %>% 
    filter(condition1 != "CSUB") %>% 
    filter(condition1 != "ASC") %>% 
    filter(Postrank != 3) %>% 
    filter(condition1 != 'ascenders') %>% 
    select(SampleID, condition1)
  
  var_info$condition1 %>%
    factor(.,levels = c("CDOM","DOM","DES")) -> group.dl
  
  dlNorm<- dlNorm[, var_info$SampleID]
  
  
  design.dl <- model.matrix(~ 0 + group.dl)
  colnames(design.dl) -> mycolnames
  
  v.dl = voom(dge.dl, design.dl, plot = F)
  vfit.dl = lmFit(v.dl, design.dl)
  
  contrast.matrix <- makeContrasts(group.dlCDOM-group.dlDOM,
                                   group.dlCDOM-group.dlDES, 
                                   group.dlDOM-group.dlDES,
                                   levels=design.dl)
  
  vfit.dl2 <- contrasts.fit(vfit.dl, contrast.matrix)
  
  efit.dl2 = eBayes(vfit.dl2)
  
  p.dl.limma2 = efit.dl2[["p.value"]]
  head(p.dl.limma2)
  
  saveRDS(v.dl, glue("manuscript/brain/manuscript70/results/limma_vdl_{my_region}_DOM2"))
  
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
    
    for(c in 1 : 3){
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
  
  saveRDS(q.dl,glue("manuscript/brain/manuscript70/results/results_RNAseqRDS/limma_eFDR_{my_region}_DOM2_cutoff5_R{R}_contrast.RDS"))
  
  q.dl <- readRDS(glue("manuscript/brain/manuscript70/results/results_RNAseqRDS/limma_eFDR_{my_region}_DOM2_cutoff5_R{R}_contrast.RDS"))
  
  
  # png(filename = glue("manuscript/brain/manuscript70/results/img/eFDR_hist_{my_region}_DOM2R{R}_contrast.png"),
  #     width = 18, height = 17, units = "cm", res = 600)
  # # hist(q.dl-p.dl.limma2, main = glue("{my_region}"))
  # invisible(dev.off())
  
  efit.dl2[["p.value"]] <- q.dl
  row.names(q.dl) <- NULL
  sum(duplicated(row.names(efit.dl2$coefficients)))
  
  tmp1 <- contrasts.fit(efit.dl2, coef = 1) # 
  tmp2 <- contrasts.fit(efit.dl2, coef = 2) # 
  tmp3 <- contrasts.fit(efit.dl2, coef = 3) # 
  
  limma_list <- list()
  
  topTable(tmp1, sort.by = "P", n = Inf) %>% 
    rownames_to_column('ensgene') %>% 
    left_join(grcm38) %>%
    filter(!is.na(symbol)) %>% 
    select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$cdom
  
  topTable(tmp2, sort.by = "P", n = Inf) %>% 
    rownames_to_column('ensgene') %>% 
    left_join(grcm38) %>%
    filter(!is.na(symbol)) %>% 
    select(symbol,logFC,P.Value,adj.P.Val)  -> limma_list$cdes
  
  
  topTable(tmp3, sort.by = "P", n = Inf) %>% 
    rownames_to_column('ensgene') %>% 
    left_join(grcm38) %>%
    filter(!is.na(symbol)) %>% 
    select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$domdes
  
  
  saveRDS(limma_list,glue("manuscript/brain/manuscript70/results/results_RNAseqRDS/limma_{my_region}_DOM2.RDS"))
  
}




limma_list<- readRDS("manuscript/brain/manuscript70/results/results_RNAseqRDS/limma_MEA_DOM.RDS") %>% 
  map(~distinct(.)) %>% 
  # map(~filter(.,abs(logFC) >= 0.2)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

#quick look at number of genes
limma_list %>% map(~filter(., P.Value<0.05)) %>% 
  map(~summarise(.,Up = sum(logFC>0.2),
                 Down = sum(logFC<0.2))) %>% 
  map(~mutate(.,Total = Up + Down))

limma_list %>% map(~hist(.$logFC))


cdom <- limma_list$cdom
cdes <- limma_list$cdes 
domdes <- limma_list$domdes 


# #mPFC
p_countdata <- read_csv("manuscript/brain/PFC_counts.csv")
colnames(p_countdata)[1] <- "ensgene"
colnames(p_countdata)[c(2:68)] <- substr(colnames(p_countdata)[c(2:68)], 7, 13)

dlNorm_list <- list(p_countdata )
names(dlNorm_list) <- c("mPFC")

#Getting metadata ready 
coldata <- read_csv("manuscript/brain/sample70min_table.csv")
head(coldata)
str(coldata)

#fixing things
coldata$region <- gsub("mPF", "PFC", coldata$region)
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
  filter(SampleName != "B1.PFCB12.3.4.trim.sam.counts") 

#just get samples I want
coldata$SampleID <- substr(coldata$SampleName, 7, 13)

coldata <- coldata %>%  select(-SampleName) %>% 
  pivot_wider(
    names_from = region, 
    values_from = region)

row.names <- coldata$SampleID
row.names(coldata) <- row.names #Assigning row names from as sample names  
head(coldata)


# dlNorm_a <- a_countdata 
# a_countdata<- a_countdata[, rownames(coldata)]
# all(rownames(coldata) == colnames(a_countdata)) #check
# 
# dlNorm_p <- p_countdata 
## p_countdata<- p_countdata[, rownames(coldata)]
# all(rownames(coldata) == colnames(p_countdata)) #check 

# dlNorm_list <- list(dlNorm_a , dlNorm_p)
# names(dlNorm_list) <- c("MeA", "mPFC")

regions <- c("mPFC")  

i = 1

# How many random sampling
R = 5000


for (i in 1:length(regions)){
  
  regions[[i]] -> my_region
  
  
  dlNorm <- dlNorm_list[[i]] %>% 
    column_to_rownames('ensgene')
  
  
  colnames(dlNorm)
  
  
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
  
  
  dge.dl$samples$group
  
  dge.dl_dom <- dge.dl[, dge.dl$samples$group %in% c("CDOM", "DOM", "DES")]
  dge.dl_dom$samples$group <- droplevels(dge.dl_dom$samples$group)
  dge.dl_dom$samples$group
  dge.dl<- dge.dl_dom
  dge.dl$samples$group
  
  var_info <- coldata %>% 
    filter(condition1 != "SUB") %>% 
    filter(condition1 != "CSUB") %>% 
    filter(condition1 != "ASC") %>% 
    filter(Postrank != 3) %>% 
    filter(condition1 != 'ascenders') %>% 
    select(SampleID, condition1)
  
  var_info$condition1 %>%
    factor(.,levels = c("CDOM","DOM","DES")) -> group.dl
  
  dlNorm<- dlNorm[, var_info$SampleID]
  
  
  design.dl <- model.matrix(~ 0 + group.dl)
  colnames(design.dl) -> mycolnames
  
  v.dl = voom(dge.dl, design.dl, plot = F)
  vfit.dl = lmFit(v.dl, design.dl)
  
  contrast.matrix <- makeContrasts(group.dlCDOM-group.dlDOM,
                                   group.dlCDOM-group.dlDES, 
                                   group.dlDOM-group.dlDES,
                                   levels=design.dl)
  
  vfit.dl2 <- contrasts.fit(vfit.dl, contrast.matrix)
  
  efit.dl2 = eBayes(vfit.dl2)
  
  p.dl.limma2 = efit.dl2[["p.value"]]
  head(p.dl.limma2)
  
  saveRDS(v.dl, glue("manuscript/brain/manuscript70/results/limma_vdl_{my_region}_DOM2"))
  
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
    
    for(c in 1 : 3){
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
  
  saveRDS(q.dl,glue("manuscript/brain/manuscript70/results/results_RNAseqRDS/limma_eFDR_{my_region}_DOM2_cutoff5_R{R}_contrast.RDS"))
  
  q.dl <- readRDS(glue("manuscript/brain/manuscript70/results/results_RNAseqRDS/limma_eFDR_{my_region}_DOM2_cutoff5_R{R}_contrast.RDS"))
  
  
  # png(filename = glue("manuscript/brain/manuscript70/results/img/eFDR_hist_{my_region}_DOM2R{R}_contrast.png"),
  #     width = 18, height = 17, units = "cm", res = 600)
  # # hist(q.dl-p.dl.limma2, main = glue("{my_region}"))
  # invisible(dev.off())
  
  efit.dl2[["p.value"]] <- q.dl
  row.names(q.dl) <- NULL
  sum(duplicated(row.names(efit.dl2$coefficients)))
  
  tmp1 <- contrasts.fit(efit.dl2, coef = 1) # 
  tmp2 <- contrasts.fit(efit.dl2, coef = 2) # 
  tmp3 <- contrasts.fit(efit.dl2, coef = 3) # 
  
  limma_list <- list()
  
  topTable(tmp1, sort.by = "P", n = Inf) %>% 
    rownames_to_column('ensgene') %>% 
    left_join(grcm38) %>%
    filter(!is.na(symbol)) %>% 
    select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$cdom
  
  topTable(tmp2, sort.by = "P", n = Inf) %>% 
    rownames_to_column('ensgene') %>% 
    left_join(grcm38) %>%
    filter(!is.na(symbol)) %>% 
    select(symbol,logFC,P.Value,adj.P.Val)  -> limma_list$cdes
  
  
  topTable(tmp3, sort.by = "P", n = Inf) %>% 
    rownames_to_column('ensgene') %>% 
    left_join(grcm38) %>%
    filter(!is.na(symbol)) %>% 
    select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$domdes
  
  
  saveRDS(limma_list,glue("manuscript/brain/manuscript70/results/results_RNAseqRDS/limma_{my_region}_DOM2.RDS"))
  
}




limma_list<- readRDS("manuscript/brain/manuscript70/results/results_RNAseqRDS/limma_mPFC_DOM2.RDS") %>% 
  map(~distinct(.)) %>% 
  # map(~filter(.,abs(logFC) >= 0.2)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

#quick look at number of genes
limma_list %>% map(~filter(., P.Value<0.05)) %>% 
  map(~summarise(.,Up = sum(logFC>0.2),
                 Down = sum(logFC<0.2))) %>% 
  map(~mutate(.,Total = Up + Down))

limma_list %>% map(~hist(.$logFC))


cdom <- limma_list$cdom
cdes <- limma_list$cdes 
domdes <- limma_list$domdes 


