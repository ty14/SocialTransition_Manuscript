#WGCNA consensus analysis 
library(limma)
library(WGCNA)
library(flashClust)
library(tidyverse)

#mPFC

dat_70 <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_MEA_ControlSUB.RDS")
datEx_70 <- t(dat_70$E)
head(datEx_70)


dat_25 <- readRDS("manuscript/brain/manuscript70/results/RDS/limma_vdl_mPFC_ControlSUB.RDS")
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
setLabels = c("mPFC_SUB", "MEA_SUB")
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
  filter(condition1 != "DES") %>%
  filter(condition1 != "DOM") %>%
  filter(condition1 != "CDOM") %>% 
  filter(region == "AMY") %>% 
  mutate(Time = "70 min") %>% 
  dplyr::select(1,3, 5:6,18:20, 7,13, 15,21) -> behav70


coldata %>% 
  filter(Postrank != 3) %>% 
  filter(Prerank != 3) %>% 
  filter(condition1 != "DES") %>%
  filter(condition1 != "DOM") %>%
  filter(condition1 != "CDOM") %>% 
  filter(region != "AMY") %>% 
  mutate(Time = "70 min") %>% 
  dplyr::select(1,3, 5:6,18:20, 7,13, 15,21) -> behav25



colnames(behav25)
colnames(behav70)
allTraits <- behav25 %>% rbind(behav70)
head(allTraits)

ifelse(allTraits$condition1 == "CSUB",1, allTraits$condition1) -> allTraits$condition1
ifelse(allTraits$condition1== "SUB", 0, allTraits$condition1) -> allTraits$condition1
ifelse(allTraits$condition1 == "ASC", -1, allTraits$condition1) -> allTraits$condition1

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
     file = "manuscript/brain/manuscript70/results/WGCNA/Consensus-dataInput_SUB.RData")

# Tutorial 2 ===========================================================
# Load the data saved in the first part
lnames = load(file = "manuscript/brain/manuscript70/results/WGCNA/Consensus-dataInput_SUB.RData")

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
my_power = 8


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
     file = "manuscript/brain/manuscript70/results/WGCNA/Consensus-NetworkConstruction-auto_power8_SUB.RData")


##Skipping Tutorial 3 for now ==========================


my_power = 8
lnames = load(file = "manuscript/brain/manuscript70/results/WGCNA/Consensus-dataInput_SUB.RData")
lnames
lnames = load(file = "manuscript/brain/manuscript70/results/WGCNA/Consensus-NetworkConstruction-auto_power8_SUB.RData")
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
lnames = load(file = "manuscript/brain/manuscript70/results/WGCNA/Consensus-dataInput_SUB.RData")
lnames
lnames = load(file = "manuscript/brain/manuscript70/results/WGCNA/Consensus-NetworkConstruction-auto_power8_SUB.RData")
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


ggplot(px2, aes(condition1, value, color = condition1))+
  geom_boxplot()+
  geom_jitter()+
  facet_wrap(~Module)

ggplot(ax2, aes(condition1, value, color = condition1))+
  geom_boxplot()+
  geom_jitter()+
  facet_wrap(~Module)

