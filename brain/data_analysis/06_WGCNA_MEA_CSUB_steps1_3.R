# libraries 
library(limma)
library(edgeR)
library(WGCNA)
library(Mus.musculus)
organism = 'org.Mm.eg.db'
library(organism, character.only = TRUE)
library(biomaRt)
library(AnnotationDbi)
library(pheatmap)
library(annotables)
library(clusterProfiler)
library(enrichplot)
library(organism, character.only = TRUE)
library(DOSE)
library(tidyverse)
grcm38 # mouse genes


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

cutoff <- 20
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
#10735   40

# Now take out groups that you want
#DOMs first 
dge.dl$samples$group

dge.dl_dom <- dge.dl[, dge.dl$samples$group %in% c("CSUB", "SUB", "ASC")]
dge.dl_dom$samples$group <- droplevels(dge.dl_dom$samples$group)
dge.dl_dom$samples$group
dge.dl<- dge.dl_dom
dge.dl$samples$group

coldata %>% 
  filter(condition1 != "DOM") %>%
  filter(condition1 != "DES") %>%
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


# Many functions expect the matrix to be transposed
datExpr <- t(v.dl$E)
## check rows/cols
nrow(datExpr)
ncol(datExpr)
rownames(datExpr)

# getting trait data
colnames(coldata)
amy <- coldata[c(3,4,9:11,14:16,19)]

ifelse(amy$condition1 == "CSUB", 1, amy$condition1) -> amy$condition2
ifelse(amy$condition2 == "SUB", 0, amy$condition2) -> amy$condition2
ifelse(amy$condition2 == "ASC", -1, amy$condition2) -> amy$condition2


amy1 <- amy %>% 
  filter(condition1 != "DES") %>%
  filter(condition1 != "DOM") %>%
  filter(condition1 != "CDOM") %>% 
  dplyr::select(-condition1)

#make everything numeric 
amy1[1:9] <- lapply(amy1[1:9], as.numeric)


sampleTree = hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)

par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 60, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 60, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExprx = datExpr[keepSamples, ]
nGenes = ncol(datExprx)
nSamples = nrow(datExprx)


 # remove B9.3.1 (ASC)

# Re-cluster samples with out any samples filtered 
sampleTree2 = hclust(dist(datExprx), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(amy1, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(amy1),
                    main = "AMY dendrogram and trait heatmap(all samples)")

Samples = rownames(datExpr)
collectGarbage()
saveRDS(datExprx,"brain/results/WGCNA/datExpr_MEA_CSUB_outliersRemoved.RDS")
saveRDS(datExpr,"brain/results/WGCNA/datExpr_MEA_CSUB.RDS")
# ===========================================================================

# Run WGCNA for each dom, descender, ascender, subordinate

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=50, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExprx, powerVector = powers, verbose = 5, )
# Plot the results:
print(sft)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;


# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main ="MEA 70 min:: Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");


# this line corresponds to using an R^2 cut-off of h

abline(h=0.85,col="red")  #changed to 0.8 (originally 0.9 in the tutorial)
abline(h=0.90,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = "MEA 70 min:Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
print(sft)


# TOMType = "signed", networkType = "signed",
# TOMType = "unsigned", networkType = "unsigned",
# TOMType = c("signed", networkType = "signed hybrid")
cor <- WGCNA::cor
WGCNA_get_net <- function(my_tissue = "MEA70min",
                          my_power =  9, 
                          my_TOMType ="signed", 
                          my_networkType = "signed hybrid"){
  
  x <- readRDS("brain/results/WGCNA/datExpr_MEA_CSUB_outliersRemoved.RDS") 
  
  net = blockwiseModules(x, 
                         power = my_power,
                         TOMType = my_TOMType, 
                         networkType = my_networkType,
                         minModuleSize = 30,
                         reassignThreshold = 0, 
                         mergeCutHeight = 0.25,
                         numericLabels = TRUE, 
                         pamRespectsDendro = FALSE,
                         saveTOMs = FALSE,
                         verbose = 3)
  
  
  saveRDS(net, "brain/results/WGCNA/net_MEA_CSUB_Power9_outliersRemoved.RDS")
  
}

# WGCNA_get_net("MEA70min", 8, "signed", "signed hybrid")
# WGCNA_get_net("MEA70min", 7, "signed", "signed hybrid")
WGCNA_get_net("MEA70min", 9, "signed", "signed hybrid")
# ========================================================================================

net <- readRDS(glue::glue("brain/results/WGCNA/net_MEA_CSUB_Power9_outliersRemoved.RDS"))

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath

dev.off() # make sure you do this before AND after 
png(file = "brain/results/img/cluster_dendo_MEA_CSUB_Power9_outliersRemoved.png",
    width=600, height=350)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "MEA CSUB Power 9")

dev.off()





# moduleLabels = net$colors
# moduleColors = labels2colors(net$colors)
MEs = net$MEs
# geneTree = net$dendrograms[[1]]

############ STEP 2


datExpr <- readRDS("brain/results/WGCNA/datExpr_MEA_CSUB_outliersRemoved.RDS") 
net <- readRDS("brain/results/WGCNA/net_MEA_CSUB_Power8_outliersRemoved.RDS")
dim(datExpr)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

moduleNumber = length(unique(moduleColors))
rownames(datExpr)


# getting trait data
coldata <- coldata %>% 
  filter(SampleID != "B9.3.1.") 


colnames(coldata)

amy <- coldata[c(3,4,9:11,14:16,18,19)]

ifelse(amy$condition1 == "CSUB", 1, amy$condition1) -> amy$condition2
ifelse(amy$condition2 == "SUB", 0, amy$condition2) -> amy$condition2
ifelse(amy$condition2 == "ASC", -1, amy$condition2) -> amy$condition2


amy1 <- amy %>% 
  filter(condition1 != "DES") %>%
  filter(condition1 != "DOM") %>%
  filter(condition1 != "CDOM") %>% 
  dplyr::select(-condition1)

#make everything numeric 
amy1[1:10] <- lapply(amy1[1:10], as.numeric)


head(amy1)

# Re-cluster samples with out any samples filtered 
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(amy1, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(amy1),
                    main = "AMY dendrogram and trait heatmap CDOM")

datTraits <- amy1



#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

saveRDS(MEs,"brain/results/WGCNA/MEA_CSUB_MEs_outliersRemoved.RDS")

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================
# # 
# sizeGrWindow(10,6)

dev.off() # make sure you do this before AND after

png(file = "brain/results/img/module_trait_MEA_CSUB_outlierRemoved.png",
    width=800, height=890, res = 130)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = "CSUB MEA 70min: Module-trait relationships")
#
#
dev.off()
# remember - they are correlation, not regression

#=====================================================================================
#
#  Code chunk - Won's code to look at regression 
#
#=====================================================================================
# Module of interest via significance 
#SUB: red, purple
#ASC: none
#CSUB: cyan
#CORT: cyan, purple, royalblue, red


ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleID") %>%
  pivot_longer(cols = 2:21, names_to = "Module") %>% 
  full_join(coldata)

head(ME_df)

ME_df$condition1
ME_df <- ME_df %>%
  mutate(status = condition1) %>% 
  filter(condition1 != "DOM")%>%
  filter(condition1 != "DES") %>%
  filter(condition1 != "CDOM")
ME_df$Module <- gsub("ME", "", ME_df$Module)

ME_df$status <- factor(ME_df$status, levels = c("CSUB", "SUB", "ASC"))

SUB1 <- ME_df %>% filter(Module == "red")
SUB2 <- ME_df %>% filter(Module == "purple")

SUB <- SUB1 %>% rbind(SUB2)


CSUB <- ME_df %>% filter(Module == "cyan")

CORT1 <- ME_df %>% filter(Module == "royalblue")

CORT <- SUB1 %>% rbind(SUB2, CORT1)


ME_df %>%
  ggplot(aes(post_Ncort, value, color = status, fill = status))+
  geom_point(size = 2, shape = 21, alpha = 0.3)+
  geom_smooth(method = "lm", se =F, alpha = 0.2, size = 1.2)+
  scale_color_manual(values = viridis::viridis(3))+
  scale_fill_manual(values = viridis::viridis(3))+
  facet_wrap(~Module, scales = "free_y")+
  theme(legend.position = "top")+
  labs(x = "CORT 70 min after reorganization",
       y = "Module eigengene") + theme_classic() -> p1
p1

ggsave(filename = "brain/results/img/MEA_eigengene_CORT_CSUB_outlierRemoved.png",
       p1,
       height = 10, width = 10, dpi = 130)

CORT %>%
  ggplot(aes(post_Ncort, value, color = status, fill = status))+
  geom_point(size = 2, shape = 21, alpha = 0.3)+
  geom_smooth(method = "lm", se =F, alpha = 0.2, size = 1.2)+
  scale_color_manual(values = viridis::viridis(3))+
  scale_fill_manual(values = viridis::viridis(3))+
  facet_wrap(~Module, scales = "free_y")+
  theme(legend.position = "top")+
  labs(x = "CORT 70 min after reorganization",
       y = "Module eigengene") + theme_classic() -> p1x
p1x

ggsave(filename = "brain/results/img/MEA_eigengene_CORT_CSUB_SIGN.png",
       p1x,
       height = 2.5, width = 6, dpi = 130)


source("functions/geom_boxjitter.R")

p2 <- ME_df %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(3))+
  scale_fill_manual(values = viridis::viridis(3))+
  facet_wrap(~Module, scales = "free_y")+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "MeA: CSUB")+ theme_classic()+ theme(legend.position = "none")

p2

ggsave(filename = "brain/results/img/MEA_CSUB_eigengene_boxplot_outliersRemoved.png",
       p2,
       height = 15, width = 15, dpi = 130)



p2x <- SUB %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(3))+
  scale_fill_manual(values = viridis::viridis(3))+
  facet_wrap(~Module, scales = "free_y")+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "MeA: SUB")+ theme_classic()+ theme(legend.position = "none")

p2x

ggsave(filename = "brain/results/img/MEA_SUB_eigengene_boxplot_SIGN.png",
       p2x,
       height = 3.5, width = 6.5, dpi = 130)


p2xx <- CSUB %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(3))+
  scale_fill_manual(values = viridis::viridis(3))+
  facet_wrap(~Module, scales = "free_y")+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "MeA: CSUB")+ theme_classic()+ theme(legend.position = "none")

p2xx
ggsave(filename = "brain/results/img/MEA_CSUB_eigengene_boxplot_SIGN.png",
       p2xx,
       height = 2.5, width = 2.5, dpi = 130)





ME_df %>%
  ggplot(aes(AggRec70min, value, color = status, fill = status))+
  geom_point(size = 2, shape = 21, alpha = 0.3)+
  geom_smooth(method = "lm", se = F,  alpha = 0.2, size = 1.2)+
  scale_color_manual(values = viridis::viridis(3))+
  scale_fill_manual(values = viridis::viridis(3))+
  facet_wrap(~Module, scales = "free_y")+
  theme(legend.position = "top")+
  labs(x = "Aggression Received after reorganization",
       y = "Module eigengene") + theme_classic() -> p3
p3

# ggsave(filename = "manuscript/brain/results/results_figures/MEA_eigengene_AggRec_ASUB.png",
# p3,
# height = 10, width = 10, dpi = 150)





#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Define variable David's score containing the David's score column of datTrait
status = as.data.frame(datTraits$condition2);
names(status) = "trait"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, status, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(status), sep="");
names(GSPvalue) = paste("p.GS.", names(status), sep="");


geneModuleMembership %>% 
  rownames_to_column("ensgene") %>% 
  left_join(geneTraitSignificance %>% rownames_to_column("ensgene")) -> gene_MM_TS_all



# to make holistic freaking dataframe ==================================================
gene_MM_TS_all$module <- moduleColors


get_mm <- function(x){
  x$moduleMembership <- x[colnames(x) == paste("MM",x$module,sep = "")] %>%
    unlist %>%
    as.numeric
  xx <- x %>%
    dplyr::select(ensgene,module,moduleMembership,GS.trait)
  return(x)
}


wgcna_whole <- get_mm(gene_MM_TS_all[1,])

for(i in 2:nrow(gene_MM_TS_all)){
  wgcna_whole <- rbind(wgcna_whole,get_mm(gene_MM_TS_all[i,]))
}

wgcna_whole %>%
  rename(GS.status = GS.trait) -> wgcna_whole


wgcna_whole %>%
  left_join(geneTraitSignificance %>% rownames_to_column("ensgene")) -> wgcna_all

saveRDS(wgcna_all, "brain/results/WGCNA/MEA_CSUB_WGCNA_MM_GS_all_outliersRemoved.RDS")




#######################################
# linear model
datExpr <- readRDS("brain/results/WGCNA/datExpr_MEA_CSUB_outliersRemoved.RDS") 
net <- readRDS("brain/results/WGCNA/net_MEA_CSUB_Power8_outliersRemoved.RDS")
MEs <- readRDS("brain/results/WGCNA/MEA_CSUB_MEs_outliersRemoved.RDS")

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes




colnames(coldata)
coldata$condition1 <- factor(coldata$condition1, levels= c("CSUB", "SUB", "ASC"))


orderMEs(MEs0) %>% 
  rownames_to_column("SampleID") %>%
  left_join(coldata) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  mutate(conditionx = factor(condition1, levels = c("SUB", "CSUB", "ASC"))) %>% 
  relocate(condition1,conditionx, batch, post_Ncort, Postds, Preds, AggGiven70min, AggRec70min ) %>% 
  dplyr::select(-SampleID, -post_idbatch, -mean_con_ng_ul, -pre_idbatch, -pre_idbatchcage,
                -time, -Prerank, -Postrank, -condition, -groupEX, -wt_d4, -wt_d8, -wt_12)-> ME_df

lm_result_list <- list()

library(lme4)
library(lmerTest)
# library(brms)
# library(tidybayes)


for(x in 1:length(MEs0)){
  k = x + 8
  ME_df[,c(1:8,k)] -> df
  md <- gsub("ME","",colnames(df)[9])
  colnames(df)[9] <- "module"
  lmer(module ~ condition1 +(1|batch) , data = df) -> mod1
  lmer(module ~ conditionx +(1|batch) , data = df) -> mod2
  summary(mod1)
  summary(mod2)
  # brm(module ~ status +(1|cohort) , data = df) -> bmod1
  # brm(module ~ statusx +(1|cohort) , data = df) -> bmod2
  # # summary(bmod1)
  # # plot(bmod1)
  # bmod1 %>% 
  #   gather_draws(b_statusSubdominant, b_statusSubordinate) %>% 
  #   median_qi() -> sum1
  # 
  # bmod2 %>% 
  #   gather_draws(b_statusxSubordinate) %>% 
  #   median_qi() -> sum2
  # 
  # rbind(sum1, sum2) %>% 
  #   as.data.frame() %>% 
  #   cbind(key = c("Alpha-subdom","Alpha-Sub","Subdom-Sub")) %>% 
  #   mutate(module = md) -> lm_result_list[[x]] 
  
  
  df
  rbind(summary(lm(df[,9] ~ df[,1]))$coefficients[2,],
        summary(lm(df[,9] ~ df[,1]))$coefficients[3,],
        summary(lm(df[,9] ~ df[,2]))$coefficients[3,])%>%
    as.data.frame() %>%
    cbind(key = c("CSUB-SUB","CSUB-ASC","SUB-ASC")) %>%
    mutate(module = md) -> lm_result_list[[x]]
}

lm_result_list %>% 
  do.call(rbind,.)%>% 
  filter(.,module != "grey") -> lm_result_all


saveRDS(lm_result_all,"brain/results/WGCNA/MEA_CSUB_lm_result_all_outliersRemoved.RDS")


#module for cort 

for(x in 1:length(MEs0)){
  k = x + 8
  ME_df[,c(1:8,k)] -> df
  md <- gsub("ME","",colnames(df)[9])
  colnames(df)[9] <- "module"
  
  lmer(module ~ post_Ncort +(1|batch) , data = df) -> mod1
  
  summary(mod1)
  
  df
  rbind(summary(lm(df[,9] ~ df[,1]))$coefficients[2,])%>%
    as.data.frame() %>%
    cbind(key = c("CORT")) %>%
    mutate(module = md) -> lm_result_list[[x]]
}

lm_result_list %>% 
  do.call(rbind,.)%>% 
  filter(.,module != "grey") -> lm_result_all

saveRDS(lm_result_all,"brain/results/WGCNA/MEA_CSUB_lm_result_all_outliersRemoved_CORT.RDS")

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================
# Module of interest via significance 
#SUB: red, purple
#ASC: none
#CSUB: cyan
#CORT: cyan, purple, royalblue, red


my_trait = "Status"
module = "red"
module = "purple"
module = "cyan"
module = 'royalblue'
module_list = c("red", "purple", "cyan", "royalblue" )
hub_gene_list = vector('list', length = length(module_list))
names(hub_gene_list) <- module_list


column = match(module, modNames);
moduleGenes = moduleColors==module;
colnames(datExpr)[moduleColors==module] -> module_gene
# 
# 
# sizeGrWindow(7, 7);
# par(mfrow = c(1,1));
# verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                    abs(geneTraitSignificance[moduleGenes, 1]),
#                    xlab = paste("Module Membership in", module, "module"),
#                    ylab = glue::glue("Gene significance for {my_trait}"),
#                    main = paste("Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
# datExpr[moduleGenes,]

my_MM_threshold = 0.8
my_GS_threshold = 0.2


for (module in module_list){
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  
  colnames(datExpr)[moduleColors==module] -> module_gene
  
  
  grcm38 %>% 
    dplyr::select(ensgene, symbol, chr, description) %>% 
    filter(ensgene %in% module_gene) -> module_gene_info
  
  
  
  geneModuleMembership %>% 
    rownames_to_column("ensgene") %>% 
    left_join(geneTraitSignificance %>% rownames_to_column("ensgene")) -> gene_MM_TS
  
  
  
  gene_MM_TS %>% 
    filter(ensgene %in% module_gene) %>% 
    dplyr::select(ensgene, glue::glue("MM{module}"), GS.trait) %>% 
    filter(abs(GS.trait) >= my_GS_threshold) -> x
  
  x[x[,glue::glue("MM{module}")]>my_MM_threshold,] -> hub_genes
  
  hub_genes %>% 
    left_join(module_gene_info) %>% 
    mutate(moduleName = glue::glue("{module}")) %>% 
    rename(moduleMembership = glue::glue("MM{module}")) -> hub_genes
  hub_gene_list[[module]] <- hub_genes
  
}



hub_gene_list %>% 
  do.call(rbind,.)%>% 
  unique() -> hubgenes_df

hubgenes_df %>% 
  dplyr::arrange(desc(GS.trait)) %>% 
  dplyr::select(-ensgene, -chr) %>% 
  head(15)

write.csv(hubgenes_df, "brain/results/WGCNA/WGCNA_tables/MEA_CSUB_hubgene_list_outliersRemoved.csv", row.names = F)






ADJ=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ, moduleColors) 

Alldegrees1 %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38 %>% dplyr::select(ensgene,symbol)) -> kIN_df

xx <- cbind(ensgene = colnames(datExpr), module = moduleColors) %>% 
  as.tibble()


kIN_df %>% 
  left_join(xx) -> kIN_df

saveRDS(kIN_df,"brain/results/WGCNA/kIN/MEA_CSUB_kIN_dataframe_outliersRemoved.RDS")

####################################################
# GO analysis for wgcna
# run step 2 first! 

# datExpr <- readRDS("manuscript/brain/manuscript70/results/wgcna/datExpr_MEA_CDOM.RDS") 
# net <- readRDS("manuscript/brain/manuscript70/results/wgcna/net_MEA_CDOM_Power9.RDS")

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

moduleNumber = length(unique(moduleColors))

modNames = substring(names(MEs), 3)


## ==============================================
gettop10GO_WGCNA <- function(module,my_showCategory = 10){
  
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  colnames(datExpr)[moduleColors==module] -> module_gene
  
  
  grcm38 %>% 
    filter(ensgene %in% module_gene) %>% 
    filter(!is.na(entrez)) %>% 
    dplyr::select(entrez) -> go_df_wgcna
  
  
  ggo <- enrichGO(gene = go_df_wgcna$entrez %>% unique(),
                  OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                  keyType = "ENTREZID",
                  ont = 'BP',
                  readable = T,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.50)
  
  
  fortify(
    ggo,
    showCategory = my_showCategory,
    by = "Count",
    split = NULL,
    includeAll = TRUE
  ) %>% 
    dplyr::arrange(desc(GeneRatio)) %>% 
    mutate(module = module) -> temp1
  
  return(rbind(temp1))
  
}

my_ont = "BP"
my_showCategory = 100


module_list %>% unique() -> allcolors

WGCNA_GOs <- vector('list', length(allcolors))

for(i in 1:length(allcolors)){
  gettop10GO_WGCNA(allcolors[i],my_showCategory) -> WGCNA_GOs[[i]]
}

WGCNA_GOs %>% 
  do.call(rbind,.) -> wgcna_all_gos

write.csv(wgcna_all_gos, 
          "brain/results/wgcna/wgcna_tables/CSUB_MEA_wgcna_all_gos_catogeryBP_outliersRemoved.csv",
          row.names = F)

# saveRDS(wgcna_all_gos,"manuscript/brain/manuscript/results/wgcna/CDOM_MEA_wgcna_all_gos.RDS")


#COME BACK LATER TO DO THE REST OF THIS IF NEEDED