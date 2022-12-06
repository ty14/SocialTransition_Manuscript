#WGCNA consensus analysis 
library(limma)
library(WGCNA)
library(flashClust)
library(tidyverse)

#mPFC

dat_70 <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_MEA_ControlDD.RDS")
datEx_70 <- t(dat_70$E)
head(datEx_70)


dat_25 <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_mPFC_ControlDD.RDS")
datEx_25 <- t(dat_25$E)
head(datEx_25)


x <- colnames(datEx_70)
y <- colnames(datEx_25)
x[x %in% y] -> both


datEx_25[,both] -> dat_25
datEx_70[,both] -> dat_70

# Tutorial 1 ===========================================================

# We work with two sets:
nSets = 2
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("mPFC_DD", "MEA_DD")
shortLabels = c("mPFC", "MEA")
# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = as.data.frame(dat_25))
names(multiExpr[[1]]$data) = colnames(dat_25)
rownames(multiExpr[[1]]$data) = rownames(dat_25)
multiExpr[[2]] = list(data = as.data.frame(dat_70))
names(multiExpr[[2]]$data) = colnames(dat_70)
rownames(multiExpr[[2]]$data) = rownames(dat_70)
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)

exprSize

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3)
gsg$allOK

if (!gsg$allOK){
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {if (sum(!gsg$goodSamples[[set]]))
    printFlush(paste("In set", setLabels[set], "removing samples",
                     paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes]
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}

sampleTrees = list()
for (set in 1:nSets){
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")}

par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7)

# bring in behav data 

coldata <- read_csv("manuscript/brain/sample70min_table.csv")
head(coldata)

coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "DOM", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "SUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "DES", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "ASC", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 4, "CSUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 1, "CDOM", coldata$condition1)


coldata %>% 
  filter(Postrank != 3) %>% 
  filter(Prerank != 3) %>% 
  filter(condition1 != "ASC") %>%
  filter(condition1 != "SUB") %>%
  filter(condition1 != "CSUB") %>% 
  filter(region == "AMY") %>% 
  mutate(Time = "70 min") %>% 
  dplyr::select(1,3, 5:6,18:20, 7,13, 15,21) -> behav70


coldata %>% 
  filter(Postrank != 3) %>% 
  filter(Prerank != 3) %>% 
  filter(condition1 != "ASC") %>%
  filter(condition1 != "SUB") %>%
  filter(condition1 != "CSUB") %>% 
  filter(region != "AMY") %>% 
  mutate(Time = "70 min") %>% 
  dplyr::select(1,3, 5:6,18:20, 7,13, 15,21) -> behav25



colnames(behav25)
colnames(behav70)
allTraits <- behav25 %>% rbind(behav70)
head(allTraits)

ifelse(allTraits$condition1 == "CDOM", 1, allTraits$condition1) -> allTraits$condition1
ifelse(allTraits$condition1== "DOM", 0, allTraits$condition1) -> allTraits$condition1
ifelse(allTraits$condition1 == "DES", -1, allTraits$condition1) -> allTraits$condition1

# Form a multi-set structure that will hold the clinical traits.
Traits = vector(mode="list", length = nSets)
for (set in 1:nSets){
  setSamples = rownames(multiExpr[[set]]$data)
  traitRows = match(setSamples, allTraits$SampleName)
  Traits[[set]] = list(data = allTraits %>% 
                         filter(SampleName %in% setSamples))
  Traits[[set]]$data <- Traits[[set]]$data %>% 
    column_to_rownames('SampleName')
}
collectGarbage()
# Define data set dimensions
nGenes = exprSize$nGenes
nSamples = exprSize$nSamples


save(multiExpr, Traits, nGenes, nSamples, setLabels, shortLabels, exprSize,
     file = "manuscript/brain/manuscript70/results/WGCNA/Consensus-dataInput_DD.RData")

# Tutorial 2 ===========================================================
# Load the data saved in the first part
lnames = load(file = "manuscript/brain/results_WGCNA/Consensus-dataInput_DD.RData")

# The variable lnames contains the names of loaded variables.
lnames
# Get the number of sets in the multiExpr structure.
nSets = checkSets(multiExpr)$nSets

# Choose a set of soft-thresholding powers
powers = c(seq(3,10,by=1), seq(12,20, by=2))
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets)
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]])
collectGarbage()
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity")
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4)
for (set in 1:nSets){
  for (col in 1:length(plotCols)){
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)

par(mfcol = c(2,2))
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7
for (col in 1:length(plotCols)) for (set in 1:nSets){
  if (set==1){
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col])
    addGrid()
  }
  if (col==1){
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set])
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set])
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) 
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) 
}

# dev.off()


# saved graphs on slack - ask becca about this
my_power = 6


net = blockwiseConsensusModules(
  multiExpr, power = my_power, minModuleSize = 50, deepSplit = 2,
  # maxBlockSize = 16000,
  pamRespectsDendro = FALSE,
  mergeCutHeight = 0.25, numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = TRUE, verbose = 5)


consMEs = net$multiMEs
moduleLabels = net$colors
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]]

sizeGrWindow(8,6)
# pdf(file = "manuscript/brain/Plots/ConsensusDendrogram-auto.pdf", wi = 8, he = 6)
plotDendroAndColors(consTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
# dev.off()




save(consMEs, moduleLabels, moduleColors, consTree,
     file = "manuscript/brain/manuscript70/results/WGCNA/Consensus-NetworkConstruction-auto_power6_DD.RData")


##Skipping Tutorial 3 for now ==========================


my_power = 6
lnames = load(file = "manuscript/brain/manuscript70/results/WGCNA/Consensus-dataInput_DD.RData")
lnames
lnames = load(file = "manuscript/brain/manuscript70/results/WGCNA/Consensus-NetworkConstruction-auto_power6_DD.RData")
lnames




exprSize = checkSets(multiExpr)
nSets = exprSize$nSets
# Set up variables to contain the module-trait correlations
moduleTraitCor = list()
moduleTraitPvalue = list()
# Calculate the correlations
for (set in 1:nSets){
  moduleTraitCor[[set]] = cor(consMEs[[set]]$data, 
                              Traits[[set]]$data, 
                              use = "p")
  moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set])
}


# Convert numerical lables to colors for labeling of modules in the plot
MEColors = labels2colors(as.numeric(substring(names(consMEs[[1]]$data), 3)))
MEColorNames = paste("ME", MEColors, sep="")
# Open a suitably sized window (the user should change the window size if necessary)
sizeGrWindow(10,7)
#pdf(file = "Plots/ModuleTraitRelationships-female.pdf", wi = 10, he = 7)
# Plot the module-trait relationship table for set number 1
set = 1
textMatrix = paste(signif(moduleTraitCor[[set]], 2), "\n(",
                   signif(moduleTraitPvalue[[set]], 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor[[set]])
par(mar = c(6, 8.8, 3, 2.2))
labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = "Consensus module--trait relationships ")
# dev.off()



# Initialize matrices to hold the consensus correlation and p-value
consensusCor = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))
consensusPvalue = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))
# Find consensus negative correlations
negative = moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0
consensusCor[negative] = pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative])
consensusPvalue[negative] = pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative])
# Find consensus positive correlations
positive = moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0
consensusCor[positive] = pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive])
consensusPvalue[positive] = pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive])

textMatrix = paste(signif(consensusCor, 2), "\n(",
                   signif(consensusPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor[[set]])
sizeGrWindow(10,7)
png(file = "ModuleTraitRelationships-consensus.png")
par(mar = c(6, 8.8, 3, 2.2))
labeledHeatmap(Matrix = consensusCor,
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Consensus module--trait relationships across\n",
                            paste(setLabels, collapse = " and ")))

dev.off()

##cyan and light green
# overlap before mPFC and MEA



# Tutorial 5 ===========================================================




nSets = 2
# For easier labeling of plots, create a vector holding descriptive names of the two sets.

my_power = 6
lnames = load(file = "manuscript/brain/manuscript70/results/WGCNA/Consensus-dataInput_DD.RData")
lnames
lnames = load(file = "manuscript/brain/manuscript70/results/WGCNA/Consensus-NetworkConstruction-auto_power6_DD.RData")
lnames

# Create a variable ds that will hold just the ds of mice in both sets
condition1= vector(mode = "list", length = nSets)
for (set in 1:nSets){
  condition1[[set]] = list(data = as.data.frame(Traits[[set]]$data$condition1))
  names(condition1[[set]]$data) = "condition1"
}
# Recalculate consMEs to give them color names
consMEsC = multiSetMEs(multiExpr, universalColors = moduleColors)
# We add the ds trait to the eigengenes and order them by consesus hierarchical clustering:
MET = consensusOrderMEs(addTraitToMEs(consMEsC, condition1))


sizeGrWindow(8,10)
#pdf(file = "Plots/EigengeneNetworks.pdf", width= 8, height = 10)
par(cex = 0.9)
plotEigengeneNetworks(MET, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
                      zlimPreservation = c(0.5, 1), xLabelsAngle = 90)

# Intramodular analysis: identifying genes with high Gene Significance and Module Membership

colnames(allTraits)

allTraits$region <- ifelse( grepl("AMY",allTraits$SampleName), 'MEA', "mPFC")

pdf <- allTraits %>% split(.$region)


mPFC <- MET[[1]]

px <- do.call(cbind, mPFC) %>% 
  tibble::rownames_to_column(.,var = "SampleName") %>% 
pivot_longer(.,cols = c(2:17, 19:25), names_to = "Module")


MEA <- MET[[2]]
ax <- do.call(cbind, MEA) %>% 
  tibble::rownames_to_column(.,var = "SampleName") %>% 
  pivot_longer(.,cols = c(2:17, 19:25), names_to = "Module")


px2 <- px %>%  full_join(behav25)
ax2 <- ax %>%  full_join(behav70)


px2 <- px2 %>% mutate(region = "mPFC")
ax2 <- ax2 %>% mutate(region = "MEA")


x$condition1 <- factor(x$condition1, levels= c("CDOM", "DOM", "DES"))

x1 <- px2 %>% full_join(ax2) %>% filter(Module == "data.MEcyan")
x <- px2 %>% full_join(ax2) %>% filter(Module == "data.MEgreen")
x2 <- px2 %>% full_join(ax2) %>% filter(Module == "data.MEmagenta")
x3 <- px2 %>% full_join(ax2) %>% filter(Module == "data.MElightgreen")
x <- px2 %>% full_join(ax2) %>% filter(Module == "data.MEgreenyellow")


m <- rbind(x1,x2,x3)

source("functions/geom_boxjitter.R")

x
x <- x %>%  mutate(conditionx = factor(condition1, levels = c("DOM", "CDOM", 'DES')))
xa <- x %>% filter(region == 'mPFC')

lmer(value ~ condition1 +(1|batch) , data = xa) -> mod1
lmer(value ~ conditionx +(1|batch) , data = xa) -> mod2
summary(mod1)
summary(mod2)



x %>%
  ggplot(aes(condition1, value, fill = condition1, color = condition1))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF","#5DC863FF", "#403891ff", "#25848EFF"))+
  scale_fill_manual(values = c("#29AF7FFF", "#440154FF", "#2D708EFF","#5DC863FF", "#403891ff", "#25848EFF"))+
   facet_wrap(~region, scales = "free_y")+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Dominant Consensus: Magenta")+ theme_classic()+ theme(legend.position = "none")

# mPFC
# scale_color_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+
#   scale_fill_manual(values = c("#5DC863FF", "#403891ff", "#25848EFF"))+





m$condition1 <- factor(m$condition1, levels= c("CDOM", "DOM", "DES"))
m$region<- factor(m$region, levels= c("MeA", "mPFC"))

# linear model


ma <- m %>% filter(region == "MeA") %>% filter(Module == "")

m %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  mutate(conditionx = factor(condition1, levels = c("DOM", "CDOM", 'DES'))) %>% 
  relocate(condition1,conditionx, batch, Postds, Preds, AggGiven70min, AggRec70min ) %>% 
  dplyr::select(-SampleName, -mean_con_ng_ul, -region,   -wt_d4, -wt_d8, -wt_12)-> ME_df

lm_result_list <- list()



library(lme4)
library(lmerTest)
# library(brms)
# library(tidybayes)





for(x in 1:length(m)){
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
    cbind(key = c("CDOM-DOM","CDOM_DES","DOM-DES")) %>%
    mutate(module = md) -> lm_result_list[[x]]
}

lm_result_list %>% 
  do.call(rbind,.)%>% 
  filter(.,module != "grey") -> lm_result_all


saveRDS(lm_result_all,"manuscript/brain/manuscript70/results/RDS/MEA_CDOM_lm_result_all.RDS")



## Get genes

nSets = 2
# For easier labeling of plots, create a vector holding descriptive names of the two sets.

my_power = 6
lnames = load(file = "manuscript/brain/manuscript70/results/WGCNA/Consensus-dataInput_DD.RData")
lnames
lnames = load(file = "manuscript/brain/manuscript70/results/WGCNA/Consensus-NetworkConstruction-auto_power6_DD.RData")
lnames

probes = names(multiExpr[[1]]$data)

consMEs.unord = multiSetMEs(multiExpr, universalColors = moduleLabels, excludeGrey = TRUE)
GS = list();
kME = list();
for (set in 1:nSets)
{
  GS[[set]] = corAndPvalue(multiExpr[[set]]$data, Traits[[set]]$data);
  kME[[set]] = corAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data);
}


GS.metaZ = (GS[[1]]$Z + GS[[2]]$Z)/sqrt(2);
kME.metaZ = (kME[[1]]$Z + kME[[2]]$Z)/sqrt(2);
GS.metaP = 2*pnorm(abs(GS.metaZ), lower.tail = FALSE);
kME.metaP = 2*pnorm(abs(kME.metaZ), lower.tail = FALSE);



GSmat = rbind(GS[[1]]$cor, GS[[2]]$cor, GS[[1]]$p, GS[[2]]$p, GS.metaZ, GS.metaP);
nTraits = checkSets(Traits)$nGenes
traitNames = colnames(Traits[[1]]$data)
dim(GSmat) = c(nGenes, 6*nTraits)
rownames(GSmat) = probes;
colnames(GSmat) = spaste(
  c("GS.set1.", "GS.set2.", "p.GS.set1.", "p.GS.set2.", "Z.GS.meta.", "p.GS.meta"),
  rep(traitNames, rep(6, nTraits)))
# Same code for kME:
kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[1]]$p, kME[[2]]$p, kME.metaZ, kME.metaP);
MEnames = colnames(consMEs.unord[[1]]$data);
nMEs = checkSets(consMEs.unord)$nGenes
dim(kMEmat) = c(nGenes, 6*nMEs)
rownames(kMEmat) = probes;
colnames(kMEmat) = spaste(
  c("kME.set1.", "kME.set2.", "p.kME.set1.", "p.kME.set2.", "Z.kME.meta.", "p.kME.meta"),
  rep(MEnames, rep(6, nMEs)))


annto <- probes %>% as.data.frame()
annto$ensgene <- annto$.
head()
symbol <- annto %>%  left_join(., grcm38 %>% dplyr::select(ensgene, symbol)) %>% select(symbol) %>% unique(.)

symbol2 = match(probes, symbol)
info = data.frame(Probe = probes, GeneSymbol = symbol2,
                  ModuleLabel = moduleLabels,
                  ModuleColor = labels2colors(moduleLabels),
                  GSmat,
                  kMEmat);
write.csv(info, file = "manuscript/brain/manuscript70/results/WGCNA/consensusAnalysis-CombinedNetworkResults_DD.csv",
          row.names = FALSE, quote = FALSE);


# info$GeneSymbol <- info %>%left_join(., grcm38 %>% dplyr::select(ensgene, symbol)) %>% select(symbol) %>% unique(.)


head(info)

MEColorNames

rownames

# cyan 17
# magenta 7
#lightgreen 16
# Recalculate consMEs to give them color names

colnames(info)

cyan <- info %>% select(Probe,p.kME.set1.ME7,Z.kME.meta.ME7, p.kME.metaME7 ) %>% 
  filter(p.kME.set1.ME7 <0.05) %>% rownames_to_column(., var = "ensgene") %>% left_join(., grcm38 %>% dplyr::select(ensgene, symbol, entrez)) %>% 
 filter(.,!is.na(entrez))  


ggo_wgcna <- enrichGO(gene = cyan$entrez ,
                      OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                      keyType = "ENTREZID",
                      ont = 'BP',
                      readable = T,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.10)
ggo_wgcna$Description[1:10]
top10go1$Description[6:10]
top10go1$Description[11:15]
top10go1$Description[16:20]
top10go1$Description[21:25]
