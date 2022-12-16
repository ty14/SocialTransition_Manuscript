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

ifelse(amy$condition1 == "CDOM", 1, amy$condition1) -> amy$condition2
ifelse(amy$condition2 == "DOM", 0, amy$condition2) -> amy$condition2
ifelse(amy$condition2 == "DES", -1, amy$condition2) -> amy$condition2


amy1 <- amy %>% 
  filter(condition1 != "ASC") %>%
  filter(condition1 != "SUB") %>%
  filter(condition1 != "CSUB") %>% 
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
abline(h = 55, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 55, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExprx = datExpr[keepSamples, ]
nGenes = ncol(datExprx)
nSamples = nrow(datExprx)


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
saveRDS(datExprx,"brain/results/WGCNA/datExpr_MEA_CDOM_outliersRemoved.RDS")
saveRDS(datExpr,"brain/results/WGCNA/datExpr_MEA_CDOM.RDS")
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
                          my_power =  6, 
                          my_TOMType ="signed", 
                          my_networkType = "signed hybrid"){
  
  x <- readRDS("brain/results/WGCNA/datExpr_MEA_CDOM_outliersRemoved.RDS") 
  
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
  
  
  saveRDS(net, "brain/results/WGCNA/net_MEA_CDOM_Power6_outliersRemoved.RDS")
  
}

# WGCNA_get_net("MEA70min", 8, "signed", "signed hybrid") #ugh probably go with 8???
# WGCNA_get_net("MEA70min", 7, "signed", "signed hybrid")
WGCNA_get_net("MEA70min", 6, "signed", "signed hybrid")
# ========================================================================================

net <- readRDS(glue::glue("brain/results/WGCNA/net_MEA_CDOM_Power6_outliersRemoved.RDS"))

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath

dev.off() # make sure you do this before AND after 
png(file = "brain/results/img/cluster_dendo_MEA_CDOM_Power6_outliersRemoved.png",
    width=600, height=350)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "MEA CDOM Power6")

dev.off()





# moduleLabels = net$colors
# moduleColors = labels2colors(net$colors)
MEs = net$MEs
# geneTree = net$dendrograms[[1]]

############ STEP 2


datExpr <- readRDS("brain/results/WGCNA/datExpr_MEA_CDOM_outliersRemoved.RDS") 
net <- readRDS("brain/results/WGCNA/net_MEA_CDOM_Power8_outliersRemoved.RDS")
dim(datExpr)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

moduleNumber = length(unique(moduleColors))
rownames(datExpr)


# getting trait data
coldata <- coldata %>% 
  filter(SampleID != "B12.1.3") %>% 
  filter(SampleID != "B11.4.4") %>% 
  filter(SampleID != "B13.1.4")


colnames(coldata)

amy <- coldata[c(3,4,9:11,14:16,18,19)]

ifelse(amy$condition1 == "CDOM", 1, amy$condition1) -> amy$condition2
ifelse(amy$condition2 == "DOM", 0, amy$condition2) -> amy$condition2
ifelse(amy$condition2 == "DES", -1, amy$condition2) -> amy$condition2


amy1 <- amy %>% 
  filter(condition1 != "ASC") %>%
  filter(condition1 != "SUB") %>%
  filter(condition1 != "CSUB") %>% 
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

saveRDS(MEs,"brain/results/WGCNA/MEA_CDOM_MEs_outliersRemoved.RDS")

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================
# # 
# sizeGrWindow(10,6)

dev.off() # make sure you do this before AND after

png(file = "brain/results/img/module_trait_MEA_CDOM_outlierRemoved.png",
    width=1000, height=1500, res = 130)

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
               main = "CDOM MEA 70min: Module-trait relationships")
#
#
dev.off()
# remember - they are correlation, not regression

#=====================================================================================
#
#  Code chunk - Won's code to look at regression 
#
#=====================================================================================
#DOM: lightsteelblue1, thriste2, blue, mediumpurple3, brown4, turquoise
#DES: pink, greenyellow, red, yellow
#CDOM: green
#CORT: lightsteelblue1, thriste2, blue, mediumpurple3, brown4, turquoise,
# black, sienna3, midnightblue, paleturquoise


ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleID") %>%
  pivot_longer(cols = 2:51, names_to = "Module") %>% 
  full_join(coldata)

head(ME_df)

ME_df$condition1
ME_df <- ME_df %>%
  mutate(status = condition1) %>% 
  filter(condition1 != "ASC") %>%
  filter(condition1 != "SUB") %>%
  filter(condition1 != "CSUB")
ME_df$Module <- gsub("ME", "", ME_df$Module)

ME_df$status <- factor(ME_df$status, levels = c("CDOM", "DOM", "DES"))


DOM1 <- ME_df %>% filter(Module == "lightsteelblue1")
DOM2 <- ME_df %>% filter(Module == "thistle2")
DOM3 <- ME_df %>% filter(Module == "blue")
DOM4 <- ME_df %>% filter(Module == "brown4")
DOM5 <- ME_df %>% filter(Module == "turquoise")
DOM6 <- ME_df %>% filter(Module == "mediumpurple3")

DOM <- DOM1 %>% rbind(DOM2, DOM3, DOM4, DOM5, DOM6)

DES1 <- ME_df %>% filter(Module == "pink")
DES2 <- ME_df %>% filter(Module == "greenyellow")
DES3 <- ME_df %>% filter(Module == "red")
DES4 <- ME_df %>% filter(Module == "yellow")

DES <- DES1 %>% rbind(DES2, DES3, DES4)

CDOM <- ME_df %>% filter(Module == "green")


CORT1 <- ME_df %>% filter(Module == "black")
CORT2 <- ME_df %>% filter(Module == "sienna3")
CORT3 <- ME_df %>% filter(Module == "midnightblue")
CORT4 <- ME_df %>% filter(Module == "paleturquoise")

CORT <- DOM1 %>% rbind(DOM2, DOM3, DOM4, DOM5, DOM6, CORT1, CORT2, CORT3,CORT4)


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

ggsave(filename = "brain/results/img/MEA_eigengene_CORT_CDOM_outlierRemoved.png",
       p1,
       height = 13, width = 15, dpi = 130)


CORT %>%
  ggplot(aes(post_Ncort, value, color = status, fill = status))+
  geom_point(size = 2, shape = 21, alpha = 0.3)+
  geom_smooth(method = "lm", se =F, alpha = 0.2, size = 1.2)+
  scale_color_manual(values = viridis::viridis(3))+
  scale_fill_manual(values = viridis::viridis(3))+
  facet_wrap(~Module, scales = "free_y", ncol =5)+
  theme(legend.position = "top")+
  labs(x = "CORT 70 min after reorganization",
       y = "Module eigengene") + theme_classic() -> p1x
p1x

ggsave(filename = "brain/results/img/MEA_eigengene_CORT_CDOM_SIGN.png",
       p1x,
       height = 5, width = 10, dpi = 130)

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
       title = "MeA: CDOM")+ theme_classic()+ theme(legend.position = "none")

p2

ggsave(filename = "brain/results/img/MEA_CDOM_eigengene_boxplot_outliersRemoved.png",
       p2,
       height = 15, width = 18, dpi = 130)




p2x <- DOM %>%
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
       title = "MeA: DOM")+ theme_classic()+ theme(legend.position = "none")

p2x

ggsave(filename = "brain/results/img/MEA_DOM_eigengene_boxplot_SIGN.png",
       p2x,
       height = 5, width = 8, dpi = 130)


p2xx <- DES %>%
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
       title = "MeA: DES")+ theme_classic()+ theme(legend.position = "none")

p2xx

ggsave(filename = "brain/results/img/MEA_DES_eigengene_boxplot_SIGN.png",
       p2xx,
       height = 5, width = 5, dpi = 130)


p2xxx <- CDOM %>%
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
       title = "MeA: CDOM")+ theme_classic()+ theme(legend.position = "none")

p2xxx

ggsave(filename = "brain/results/img/MEA_CDOM_eigengene_boxplot_SIGN.png",
       p2xxx,
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

saveRDS(wgcna_all, "brain/results/WGCNA/MEA_CDOM_WGCNA_MM_GS_all_outliersRemoved.RDS")




#######################################
# linear model
datExpr <- readRDS("brain/results/WGCNA/datExpr_MEA_CDOM_outliersRemoved.RDS") 
net <- readRDS("brain/results/WGCNA/net_MEA_CDOM_Power8_outliersRemoved.RDS")
MEs <- readRDS("brain/results/WGCNA/MEA_CDOM_MEs_outliersRemoved.RDS")

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes




colnames(coldata)
coldata$condition1 <- factor(coldata$condition1, levels= c("CDOM", "DOM", "DES"))


orderMEs(MEs0) %>% 
  rownames_to_column("SampleID") %>%
  left_join(coldata) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  mutate(conditionx = factor(condition1, levels = c("DOM", "CDOM", 'DES'))) %>% 
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
    cbind(key = c("CDOM-DOM","CDOM-DES","DOM-DES")) %>%
    mutate(module = md) -> lm_result_list[[x]]
}

lm_result_list %>% 
  do.call(rbind,.)%>% 
  filter(.,module != "grey") -> lm_result_all


saveRDS(lm_result_all,"brain/results/WGCNA/MEA_CDOM_lm_result_all_outliersRemoved.RDS")


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

saveRDS(lm_result_all,"brain/results/WGCNA/MEA_CDOM_lm_result_all_outliersRemoved_CORT.RDS")

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================
# Module of interest via significance 
#DOM: lightsteelblue1, thriste2, blue, mediumpurple3, brown4, turquoise
#DES: pink, greenyellow, red, yellow
#CDOM: green
#CORT: lightsteelblue1, thriste2, blue, mediumpurple3, brown4, turquoise,
# black, sienna3, midnightblue, paleturquoise


my_trait = "Status"
module = "red"
module = "lightsteelblue1"
module = "pink"
module = "brown4"
module = "green"
module = "yellow"
module = "greenyellow"
module = "thistle2"
module = 'black'
module = "mediumpurple3"
module = "paleturquoise"
module = 'blue'
module = "turquoise"
module = 'sienna3'
module = "midnightblue"
module_list = c("red", "lightsteelblue1", "pink", "brown4", "green", "yellow", "greenyellow", "thistle2", "black", 
                "mediumpurple3", "paleturquoise", "blue", "turquoise", "sienna3", "midnightblue")
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

write.csv(hubgenes_df, "brain/results/WGCNA/WGCNA_tables/MEA_CDOM_hubgene_list_outliersRemoved.csv", row.names = F)






ADJ=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ, moduleColors) 

Alldegrees1 %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38 %>% dplyr::select(ensgene,symbol)) -> kIN_df

xx <- cbind(ensgene = colnames(datExpr), module = moduleColors) %>% 
  as.tibble()


kIN_df %>% 
  left_join(xx) -> kIN_df

saveRDS(kIN_df,"brain/results/WGCNA/kIN/MEA_CDOM_kIN_dataframe_outliersRemoved.RDS")

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
          "brain/results/wgcna/wgcna_tables/CDOM_MEA_wgcna_all_gos_catogeryBP_outliersRemoved.csv",
          row.names = F)

# saveRDS(wgcna_all_gos,"manuscript/brain/manuscript/results/wgcna/CDOM_MEA_wgcna_all_gos.RDS")


#COME BACK LATER TO DO THE REST OF THIS IF NEEDED













#############STEP3
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 9);
# saveRDS(dissTOM,"manuscript/brain/results/wgcna/SUB_ASC/dissTom_MEA70minASUB.RDS")
# dissTOM <- readRDS("manuscript/brain/results/wgcna/SUB_ASC/dissTom_MEA70minASUB.RDS")

restGenes= (moduleColors != "grey")
dissTOM = 1-TOMsimilarityFromExpr(datExpr[,restGenes],
                                  power = my_power,
                                  networkType = "signed hybrid");
plotTOM = dissTOM^7;
# look at the network w/out "grey" modules
colnames(dissTOM) = rownames(dissTOM) = colnames(datExpr[restGenes])
hier1=flashClust(as.dist(dissTOM), method = "average" )



plotDendroAndColors(hier1, colors = data.frame(moduleColors[restGenes]),
                    c("Module Tree Cut"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
                    guideHang = 0.05, main = "Gene dendrogram and module colors")

# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
diag(dissTOM) = NA






png(glue("manuscript/brain/manuscript70/results/results_figures/MEA_CDOM_TOMplotsignedhybird_9.png"),   width = 12, height = 12, units = "cm", res = 600)

TOMplot(dissim = 1-dissTOM^7, 
        hier1, 
        as.character(moduleColors[restGenes]), 
        min = "Network heatmap plot, removed unassigned genes")
dev.off()




# lm_result_all %>% 
#   rename(Estimate = value) %>% 
#   mutate(my_alpha = ifelse((.lower)*(.upper) >0, 1, 0)) %>% 
#   mutate(my_alpha2 = Estimate ) -> heatmap_df

heatmap_df <- lm_result_all %>% 
  mutate(my_alpha = ifelse(`Pr(>|t|)` < 0.05, 1, 0)) %>% 
  mutate(my_alpha2 = Estimate )

moduleColors %>% 
  table() %>% 
  as.data.frame() %>% arrange(Freq)  -> modnum 
colnames(modnum) <- c("module","count")

heatmap_df %>% as_tibble() %>% 
  left_join(modnum) -> heatmap_dfx



heatmap_df %>% 
  filter(my_alpha > 0.95) %>% 
  .$Estimate %>% 
  abs() %>% 
  max() -> my_limit
rev(levels(as.factor(heatmap_df$module))) -> xx
xx[xx != "grey"] -> y_limit


heatmap_dfx %>% 
  as_tibble() %>% 
  select(module, count) %>% 
  distinct() %>% 
  # mutate(module = factor(module)) %>% 
  arrange(desc(count)) %>% 
  filter(module!="grey") %>% 
  .$module -> my_module_level


key_level <- c("CDOM-DOM","CDOM-DES","DOM-DES")


vv = length(my_module_level)-1
vlines = c(1:vv)+0.495
hlines = c(1:2)+0.51
ggplot(data = heatmap_df %>%  filter(my_alpha >0),
       aes(x = module, y =key, fill = Estimate))+
  geom_tile(color = "white",size =1 )+
  geom_vline(xintercept = vlines, color = "grey",size = 0.2)+
  geom_hline(yintercept = hlines, color = "grey",size = 0.2)+
  scale_x_discrete(limits = my_module_level, position = 'top')+
  scale_y_discrete(limits = rev(key_level))+
  scale_fill_distiller(palette = "PiYG",
                       limit = c(-my_limit,my_limit))+
  labs(x = "", y = "")+
  theme_bw(base_size = 15)+
  theme(legend.title = element_text(size =rel(1.2)),
        legend.text = element_text(size =rel(1.)),
        legend.position = "bottom", 
        legend.justification = c(0, 0),
        legend.direction = "horizontal",
        # legend.key.height = unit(1.5, 'cm'), 
        legend.key.width = unit(2, 'cm'),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "grey"),
        # panel.background = element_blank(),
        # panel.grid.minor = element_line(color = 'grey50'),
        panel.grid = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"),
        axis.text.y = element_text(hjust=0.95,vjust=0.2,size = rel(1.2)),
        axis.text.x = element_text(angle = 45,hjust=0.08,vjust=0.15,size = rel(1))) -> traitmodule2


traitmodule2



png(filename = "manuscript/brain/manuscript70/results/results_figures/traitmodule_MEA_DOMgroups.png",
    width = 18, height = 13, units = "cm", res = 600)

traitmodule2

invisible(dev.off())


#=============================================
modnum$module <- factor(modnum$module,levels=unique(modnum$module[order(modnum$count)]))
mycol <- modnum %>% 
  .$module %>% as.character() # for some reason this has to be character again to be ordered properly in the figure...!! 


modnum %>%
  filter(module != "grey") %>% 
  ggplot(aes(y = module, x = count, fill = module))+
  geom_bar(stat = 'identity', color = NA)+
  geom_text(aes(label = count),hjust = -0.2, size = 4)+
  xlim(0, max(modnum$count)+75)+
  scale_fill_manual(values = mycol)+
  theme_minimal(base_size = 7)+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(margin = margin(r=-6)),
        panel.grid = element_blank())+
  labs(x = "",
       y = "",title = "Module Size")+
  theme_minimal(base_size = 15)+
  theme(legend.position = "none", 
        legend.title = element_text(size =rel(1.2)),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size =rel(1.5)),
        legend.key.size = unit(.9, 'cm'),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_blank(),
        # panel.background = element_blank(),
        # panel.grid.minor = element_line(color = 'grey50'),
        panel.grid = element_blank(),
        axis.text.y = element_text(hjust=1,vjust=0.25,size = rel(1.25))) -> temp_p
temp_p

ggsave(filename = "manuscript/brain/manuscript70/results/results_figures/modulesize_MEA_CDOM.png",temp_p,
       width = 30, height = 14.5, units = "cm", dpi = 150)



# GO over-representation analysis ================================================
dev.off()

ggo_wgcna <- enrichGO(gene = go_df_wgcna$entrez ,
                      OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                      keyType = "ENTREZID",
                      ont = 'BP',
                      readable = T,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.10)



png(glue::glue("manuscript/brain/results/results_figures/enrichGO_MEA_70minASUB_WGCNA_{module}_BP.png"),
    width = 9, height = 4, units = 'in', res = 300)

clusterProfiler::dotplot(ggo_wgcna, orderBy = "Count")+
  ggtitle(glue::glue("MEA 70min DOM:DES: {module} Module"))+
  theme(plot.title = element_text(size = 12))

dev.off()

# KEGG over-representation enrichment analysis ================================================

kk_wgcna <- enrichKEGG(gene    = go_df_wgcna$entrez,
                       organism     = 'mmu',
                       pvalueCutoff = 0.05)
head(kk_wgcna)

png(glue("results_figures/gseKEGG_{my_tissue}_WGCNA_{module}.png"), 
    width = 9, height = 4, units = 'in', res = 300)
clusterProfiler::dotplot(kk_wgcna)+
  ggtitle(glue("KEGG - {my_tissue}: Module {module}"))+
  theme(plot.title = element_text(size = 12))+
  scale_color_viridis()

dev.off()

# entrez list for genewalk =========================================


moduleColors %>% unique() -> allcolors
for(i in 1:length(allcolors)){
  allcolors[i] -> module
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  colnames(datExpr)[moduleColors==module] -> module_gene
  
  
  grcm38 %>% 
    filter(ensgene %in% module_gene) %>% 
    filter(!is.na(entrez)) %>% 
    dplyr::select(entrez) %>% 
    write.table(glue::glue("manuscript/brain/results/results_GeneWalk/forGeneWalk_MEA70min_ASUB_WGCNA_{module}.csv"), sep = "\t",
                quote = F,
                col.names = F,
                row.names = F)
}



