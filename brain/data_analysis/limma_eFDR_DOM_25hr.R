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
dlNorm <-  read.csv("brain/AMY_25hcounts.csv", row.names = 1)
#remove zeros
dlNorm <- dlNorm[apply(dlNorm[], 1, function(x) !all(x==0)),]
#trim sample ids
colnames(dlNorm)[c(1:24)] <- substr(colnames(dlNorm)[c(1:24)], 6, 12)


# Getting metadata ready 
coldata <- read_csv("brain/sample25hr_table.csv")
head(coldata)
str(coldata)

#fixing things
coldata$region <- gsub("A", "MeA",coldata$region)
coldata$region <- gsub("P", "mPFC",coldata$region)

#fixing postbatchID 3-3Batch1
coldata$condition[is.na(coldata$condition)] <- "same"



coldata$groupEX <- coldata$group


# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))

coldata <-coldata %>% mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))
coldata <- coldata %>%  dplyr::select(-group, -period)
coldata$post_Ncort[is.na(coldata$post_Ncort)] <-  mean(coldata$post_Ncort, na.rm = TRUE)


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
colnames(coldata)

coldata$SampleID <- substr(coldata$SampleNames, 6, 12)


coldata <- coldata %>% select(-region, -SampleNames) %>% unique(.)


row.names <- coldata$SampleID
row.names(coldata) <- row.names #Assigning row names from as sample names  
head(coldata)


#check before normalizing 
dlNorm <- dlNorm[, rownames(coldata)]
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
#15127    24

# Now take out groups that you want
#DOMs first 
dge.dl$samples$group

dge.dl_dom <- dge.dl[, dge.dl$samples$group %in% c("DOM", "DES")]
dge.dl_dom$samples$group <- droplevels(dge.dl_dom$samples$group)
dge.dl_dom$samples$group
dge.dl<- dge.dl_dom
dge.dl$samples$group

coldata %>% 
  filter(condition1 != "ASC") %>%
  filter(condition1 != "SUB") %>%
  dplyr::select(SampleID, condition1,post_Ncort) -> var_info  

row.names <- var_info$SampleID

row.names(var_info) <- row.names #Assigning row names from as sample names  
head(var_info )

dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check

##following Won's code
ifelse(var_info$condition1 == "DOM", 1,-1) -> group.dl

var_info$post_Ncort %>% 
  scale %>% 
  as.numeric -> cort.dl 

design.dl <- model.matrix(~ group.dl*cort.dl)
colnames(design.dl) -> mycolnames


v.dl = voom(dge.dl, design.dl, plot = F)
vfit.dl = lmFit(v.dl, design.dl)
efit.dl = eBayes(vfit.dl)
p.dl.limma = efit.dl[["p.value"]]

saveRDS(v.dl, glue("brain/results/RDS/25hr/limma_vdl_MeA_DOM25.RDS"))

# How many random sampling
R = 5000

p.dl.rand = vector('list',length = R)

for(i in 1 : R){
  print(paste("Starting on Permutation", i))
  
  # Randomize the traits
  cort.dl.rand = sample(cort.dl)
  group.dl.rand = sample(group.dl)
  
  # Model
  design.dl.rand = model.matrix(~group.dl.rand*cort.dl.rand)
  colnames(design.dl.rand) <- mycolnames
  
  # Calculate p-values based on randomized traits
  v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
  vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
  efit.dl.rand = eBayes(vfit.dl.rand)
  p.dl.rand[[i]] = efit.dl.rand[["p.value"]]
}

q.dl = matrix(0, nrow = nrow(p.dl.limma), ncol = ncol(p.dl.limma))

for(i in 1 : R){
  print(paste("Calculating Permutation", i))
  
  temp = p.dl.rand[[i]]
  
  for(c in 1 : 4){
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

saveRDS(q.dl,glue("brain/results/RDS/25hr/limma_eFDR_MeA_DOM25_cutoff5_2000_tworand.RDS"))


replace p value with q value 
efit.dl[["p.value"]] <- q.dl
row.names(q.dl) <- NULL
sum(duplicated(row.names(efit.dl$coefficients)))

tmp_status <- contrasts.fit(efit.dl, coef = 2) # group
tmp_cort <- contrasts.fit(efit.dl, coef = 3) # cort
tmp_interaction <- contrasts.fit(efit.dl, coef = 4) # interaction

limma_list <- list()

topTable(tmp_status, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$status


topTable(tmp_cort, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  select(symbol,logFC,P.Value,adj.P.Val)  -> limma_list$cort


topTable(tmp_interaction, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$interaction


saveRDS(limma_list,"brain/results/RDS/25hr/limma_MeA_DOM25.RDS")
