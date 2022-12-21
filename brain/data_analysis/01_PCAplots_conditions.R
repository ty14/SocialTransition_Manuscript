library(limma)
library(Mus.musculus)
library(DESeq2)
library(tidyverse)


## Getting PCA from the DEG results not LIMMA

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


#just get samples I want
coldata$SampleID <- substr(coldata$SampleName, 7, 13)

# coldata <- coldata %>%  filter(post_idbatch!= "3-4Batch12") for mPFC

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

cutoff <- 50
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
# 6235  39
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
  dplyr::select(SampleID, condition1) -> var_info  

row.names <- var_info$SampleID

row.names(var_info) <- row.names #Assigning row names from as sample names  
head(var_info )

dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check

##following Won's code
var_info$condition1 %>%
  factor(.,levels = c("CDOM","DOM","DES")) -> group.dl


design.dl <- model.matrix(~ 0 + group.dl)
colnames(design.dl) -> mycolnames

v = voom(dge.dl, design.dl, plot = F)

colData <- v$design %>% data.frame()
colData <- var_info %>% rownames_to_column(var = "SampleName") %>% cbind(colData) 

rv <- rowVars((v$E))
select <- order(rv, decreasing = TRUE)[1:100]
pca1 <- prcomp(t((v$E)[select, ]))

condition3.df <- as.data.frame(colData[, 'condition1', 
                                       drop = FALSE])




d1 <- data.frame(PC1 = pca1$x[, 1], PC2 = pca1$x[, 2], PC3 = pca1$x[, 3], PC4 = pca1$x[, 4], 
                 condition3.df, name = colnames(v$E))

d1$condition1

head(d1)

pv1 <- ((pca1$sdev^2) / (sum(pca1$sdev^2)))*100
barplot(pv1, cex.names=1, xlab=paste("Principal component (PC), 1-", length(pca1$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))

var1 <- ((pca1$sdev[1]^2) / (sum(pca1$sdev^2)))*100
var2 <- ((pca1$sdev[2]^2) / (sum(pca1$sdev^2)))*100
var3 <- ((pca1$sdev[3]^2) / (sum(pca1$sdev^2)))*100
var4 <- ((pca1$sdev[4]^2) / (sum(pca1$sdev^2)))*100

## actually plots 

xlab=paste("PC1, ", round(pv1[1], 2), "%")
ylab=paste("PC2, ", round(pv1[2], 2), "%")

xlab3=paste("PC3, ", round(pv1[3], 2), "%")
ylab4=paste("PC4, ", round(pv1[4], 2), "%")
d1$condition1

ggplot(d1, aes(PC1, PC2, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = xlab, y = ylab)+
  scale_color_manual(values = viridis::viridis(3))+
  theme_classic()

ggplot(d1, aes(PC1, PC3, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = xlab, y = xlab3)+
  scale_color_manual(values = viridis::viridis(3))+
  theme_classic()

ggplot(d1, aes(PC2, PC3, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = ylab, y = xlab3)+
  scale_color_manual(values = viridis::viridis(3))+
  theme_classic()

ggplot(d1, aes(PC3, PC4, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = xlab3, y = ylab4)+
  scale_color_manual(values = viridis::viridis(3))+
  theme_classic()

ggplot(d1, aes(PC1, PC4, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = xlab, y = ylab4)+
  scale_color_manual(values = viridis::viridis(3))+
  theme_classic()


ggplot(d1, aes(PC2, PC4, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = ylab, y = ylab4)+
  scale_color_manual(values = viridis::viridis(3))+
  theme_classic()


#########
# subs
# Expression values
dlNorm <-  read.csv("brain/PFC_counts.csv", row.names = 1)
#remove zeros - Becca told me to do this since zeros can fuck things up. 
dlNorm <- dlNorm[apply(dlNorm[], 1, function(x) !all(x==0)),]

#Lines 20 to 70 is just me getting my group variables together
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

#normalize and filter with all groups -
# I normalize with all groups since I want to compare them after

dlNorm <- dlNorm[!is.na(rowSums(dlNorm)),]

d = apply(dlNorm, 2, as.numeric)
dim(d)

d0= DGEList(d, group = coldata$condition1)
dim(d0)
rownames(d0) <- rownames(dlNorm)
d0 <- calcNormFactors(d0)

cutoff <-50
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
#6571   40

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
  dplyr::select(SampleID, condition1) -> var_info  

row.names <- var_info$SampleID

row.names(var_info) <- row.names #Assigning row names from as sample names  
head(var_info )

dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check

##following Won's code
var_info$condition1 %>%
  factor(.,levels = c("CSUB","SUB","ASC")) -> group.dl


design.dl <- model.matrix(~ 0 + group.dl)
colnames(design.dl) -> mycolnames

v.dl = voom(dge.dl, design.dl, plot = F)








colData <- v$design %>% data.frame()
colData <- var_info %>% rownames_to_column(var = "SampleName") %>% cbind(colData) 

rv <- rowVars((v$E))
select <- order(rv, decreasing = TRUE)[1:100]
pca1 <- prcomp(t((v$E)[select, ]))

condition3.df <- as.data.frame(colData[, 'condition1', 
                                       drop = FALSE])




d1 <- data.frame(PC1 = pca1$x[, 1], PC2 = pca1$x[, 2], PC3 = pca1$x[, 3], PC4 = pca1$x[, 4], 
                 condition3.df, name = colnames(v$E))

d1$condition1

head(d1)

pv1 <- ((pca1$sdev^2) / (sum(pca1$sdev^2)))*100
barplot(pv1, cex.names=1, xlab=paste("Principal component (PC), 1-", length(pca1$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))

var1 <- ((pca1$sdev[1]^2) / (sum(pca1$sdev^2)))*100
var2 <- ((pca1$sdev[2]^2) / (sum(pca1$sdev^2)))*100
var3 <- ((pca1$sdev[3]^2) / (sum(pca1$sdev^2)))*100
var4 <- ((pca1$sdev[4]^2) / (sum(pca1$sdev^2)))*100

## actually plots 

xlab=paste("PC1, ", round(pv1[1], 2), "%")
ylab=paste("PC2, ", round(pv1[2], 2), "%")

xlab3=paste("PC3, ", round(pv1[3], 2), "%")
ylab4=paste("PC4, ", round(pv1[4], 2), "%")
d1$condition1

ggplot(d1, aes(PC1, PC2, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = xlab, y = ylab)+
  scale_color_manual(values = viridis::viridis(3))+
  theme_classic()

ggplot(d1, aes(PC1, PC3, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = xlab, y = xlab3)+
  scale_color_manual(values = viridis::viridis(3))+
  theme_classic()

ggplot(d1, aes(PC2, PC3, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = ylab, y = xlab3)+
  scale_color_manual(values = viridis::viridis(3))+
  theme_classic()

ggplot(d1, aes(PC3, PC4, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = xlab3, y = ylab4)+
  scale_color_manual(values = viridis::viridis(3))+
  theme_classic()

ggplot(d1, aes(PC1, PC4, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = xlab, y = ylab4)+
  scale_color_manual(values = viridis::viridis(3))+
  theme_classic()


ggplot(d1, aes(PC2, PC4, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = ylab, y = ylab4)+
  scale_color_manual(values = viridis::viridis(3))+
  theme_classic()
