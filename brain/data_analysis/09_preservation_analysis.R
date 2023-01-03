
# Load the package
library(WGCNA)
library(tidyverse)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)


domExp <- readRDS("brain/results/wgcna/datExpr_MeA_CDOM_outliersRemoved.RDS")

domColor <-readRDS("brain/results/wgcna/kIN/MEA_CDOM_kIN_dataframe_outliersRemoved.RDS")
domColor <- domColor %>% dplyr::select(ensgene, module) %>% dplyr::rename(.,geneID = ensgene)%>% dplyr::rename(.,Module = module)
  
subExp <- readRDS("brain/results/wgcna/datExpr_MeA_CSUB_outliersRemoved.RDS")

subColor <-readRDS("brain/results/wgcna/kIN/MEA_CSUB_kIN_dataframe_outliersRemoved.RDS")
subColor <- subColor %>% dplyr::select(ensgene, module) %>% dplyr::rename(.,geneID = ensgene) %>% dplyr::rename(.,Module = module)


#match genes in dom and sub expression data 
dom2sub = match(colnames(domExp), colnames(subExp));
table(is.finite(dom2sub))
domExp = domExp[,dom2sub];
all.equal(colnames(domExp), colnames(subExp))

subColor = subColor[subColor$geneID %in% domColor$geneID,]
all.equal(domColor$geneID, subColor$geneID)



# # Calculating topological overlap and clustering
library(flashClust)
#first doms
dissTOM_dom = 1-TOMsimilarityFromExpr(domExp[,domColor$geneID], power = 5, TOMType = "unsigned")
tree_dom = flashClust(as.dist(dissTOM_dom), method = "a")


#dendrogram
plotDendroAndColors(tree_dom, colors = domColor$Module,
                    c("Module Tree Cut"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
                    guideHang = 0.05, main = "DOM: gene dendrogram and module colors ")


#then subs
dissTOM_sub = 1-TOMsimilarityFromExpr(subExp[,subColor$geneID], power = 6, TOMType = "unsigned")
tree_sub = flashClust(as.dist(dissTOM_sub), method = "a")

#both dendrogram
sizeGrWindow(11, 7)
# Set up appropriate screen sectioning
layout(matrix(c(1:4), 4, 1), heights = rep(c(0.8, 0.2), 2));
# Plot the female dendrogram
plotDendroAndColors(tree_dom,domColor$Module ,
                    "DOM modules", main = "DOM gene dendrogram and module colors",
                    dendroLabels = FALSE, setLayout = FALSE, marAll = c(2,10,3,0), cex.colorLabels = 1.4,
                     cex.main = 2, cex.lab = 1.4, cex.axis = 1.4,
                    addGuide = TRUE);
# Plot the SUB dendrogram with DOM module colors
plotDendroAndColors(tree_sub, domColor$Module,
                    "DOM modules", main = "SUB gene dendrogram and DOM module colors",
                    dendroLabels = FALSE, setLayout = FALSE, marAll = c(2,10,3,0), cex.colorLabels = 1.4,
                     cex.main = 2, cex.lab = 1.4, cex.axis = 1.4,
                    addGuide = TRUE)
dev.off()

#  Module preservation

setLabels = c("DOM", "SUB")
multiExpr = list(DOM = list(data = domExp[,domColor$geneID]), SUB = list(data = subExp[,subColor$geneID]))
multiColor = list(DOM = domColor$Module, SUB= domColor$Module)
nSets = 2

#takes 2.5 hours - make use you save out put. 
# Calculate module preservation statistics
system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = c(1:2),
                          loadPermutedStatistics = FALSE,
                          networkType = "signed hybrid",
                          nPermutations = 1000,
                          verbose = 3)
} );
# Save the results
saveRDS(mp,"brain/results/RDS/MEA_modulePreservation.RDS")

mp <- readRDS("brain/results/RDS/MEA_modulePreservation.RDS") #dom and sub colors
# get observed statistics and their Z scores

ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])
# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )


# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5)
# png("manuscript/brain/manuscript70/results/results_figures/mPFC_DOMonly-modulePreservation-Zsummary-medianRank.png", wi=10, h=5,units = 'in', res = 300)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
# If plotting into a file, close it
dev.off()

##PLOT2
# Re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1];
# Exclude improper modules
plotMods = !(modColors %in% c("grey", "gold"));
# Create numeric labels for each module
labs = match(modColors[plotMods], standardColors(50));
# Start the plot: open a suitably sized graphical window and set sectioning and margins. Alternatively,
# plot into a pdf file.
# sizeGrWindow(10, 9);
png("brain/results/img/MEA_DOMonly-PreservationZStatistics.png", wi=10, h=9,units = 'in', res = 300)
par(mfrow = c(4,5))
par(mar = c(3,3,2,1))
par(mgp = c(1.6, 0.4, 0));
for (s in 1:ncol(statsZ))
{
  min = min(statsZ[plotMods, s], na.rm = TRUE);
  max = max(statsZ[plotMods, s], na.rm = TRUE);
  if (min > -max/12) min = -max/12
  plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
       main = colnames(statsZ)[s],
       cex = 2.2,
       ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
       ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
       xlim = c(30, 1000),
       cex.lab = 1.2, cex.axis = 1.2)
  labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 1, offs = 0.06);
  # text(moduleSizes[-1], statsZ[-c(1:2), s], labels = letter[-c(1:2)], col = "black"); #modColors[-2]);
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
}
# If plotting into a file, close it, otherwise it is unreadable.
dev.off();

data.frame(color = modColors[plotMods], label = labs)




#########################
#this can take several days to run. 
# library(clusterRepro)
# nSets = 2
# doClusterRepro = TRUE;
# if (doClusterRepro)
# {
#   cr = list();
#   library(clusterRepro);
#   set.seed(20);
#   for (ref in 1:nSets)
#   {
#     eigengenes = multiSetMEs(multiExpr, universalColors = multiColor[[ref]]);
#     cr[[ref]] = list();
#     for (set in 1:nSets)
#     {
#       printFlush(paste("Working on reference set", ref, " and test set", set));
#       rownames(eigengenes[[set]]$data) = rownames(multiExpr[[set]]$data);
#       print(system.time({ cr[[ref]][[set]] = clusterRepro(Centroids = as.matrix(eigengenes[[set]]$data),
#                                                           New.data = as.matrix(multiExpr[[set]]$data),
#                                                           Number.of.permutations = 2000)}))
#     }
#   }
# }
# # Save the results so they can be re-used if necessary.
# saveRDS(cr, file = "manuscript/brain/manuscript70/results/RDS/MEA_modulePreservation-cr.RDS")
# statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
# statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])

# x <- mp$preservation$log.p$ref.DOM$inColumnsAlsoPresentIn.SUB$log.p.cor.kME
# y <- mp$preservation$log.p$ref.SUB$inColumnsAlsoPresentIn.DOM$log.p.cor.kME
# km <- data.frame(x,y)
# 
# ggplot(km,aes(x,y))+
#   geom_point()


#####################I have no idea what I am doing 
hub <- read_csv("manuscript/brain/manuscript70/results/wgcna/wgcna_table/mPFC_CDOM_hubgene_list.csv")
head(hub)

hub$moduleName

#blue
b_hub <- hub %>%  filter(moduleName == "green") %>% arrange(symbol) %>% unique(.)
b <- b_hub$symbol

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_mPFC_ControlDD.RDS") %>% 
  map(~distinct(.)) %>% 
  # map(~filter(.,abs(logFC) >= .2)) %>%
  # map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol))) %>% 
  map(~dplyr::select(.,symbol)) 

x <- limma_list %>%   map2_df(.,names(.), ~mutate(.x, contrast = .y))

x2 <- x %>% dplyr::select(.,symbol) %>% unique(.)
x2 <- t(x2)

br_sign <- x2[x2 %in% b] 
brx <- as.data.frame(br_sign) %>% dplyr::select(.,symbol = br_sign) %>% arrange(symbol)
xx <- limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_MEA_ControlDD.RDS") %>% 
  map(~distinct(.)) %>% 
  # map(~filter(.,abs(logFC) >= .2)) %>%
  # map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol))) %>% 
  map2_df(.,names(.), ~mutate(.x, contrast = .y)) %>% unique(.)


br_doc <- xx[xx$symbol %in% brx$symbol,]

b_hub <- b_hub[b_hub$symbol %in% brx$symbol,]

bb <- br_doc %>% filter(contrast == "domdes") %>% arrange(symbol) %>% select(logFC, P.Value)
bb1 <-  br_doc %>% filter(contrast == "controldom") %>% arrange(symbol)%>% select(logFC, P.Value)
bb2 <- br_doc %>% filter(contrast == "controldes") %>% arrange(symbol)%>% select(logFC, P.Value)

x <- b_hub %>% select( symbol,moduleMembership) %>% cbind(bb1,bb2,bb)
write_csv(x, "manuscript/brain/manuscript70/results/MEA_black.csv")
####### subordinate data 

hub <- readRDS("manuscript/brain/manuscript70/results/wgcna/mPFC_CSUB_WGCNA_MM_GS_all.RDS")
head(hub)

subhub <- grcm38 %>% 
  dplyr::select(ensgene, symbol, chr, description) %>% 
  filter(ensgene %in% hub$ensgene) %>% full_join(hub)


subhub2 <- subhub[subhub$symbol %in% brx$symbol,] %>% arrange(symbol) %>% unique(.)


xx <- limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_mPFC_ControlSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  # map(~filter(.,abs(logFC) >= .2)) %>%
  # map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol))) %>% 
  map2_df(.,names(.), ~mutate(.x, contrast = .y)) %>% unique(.)


doc <- xx[xx$symbol %in% subhub2$symbol,]



sub <- doc %>% filter(contrast == "controlsub") %>% arrange(symbol)%>% select(logFC, P.Value)
sub1 <- doc %>% filter(contrast == "controlsub") %>% arrange(symbol)%>% select(logFC, P.Value)
sub2 <- doc %>% filter(contrast == "controlsub") %>% arrange(symbol)%>% select(logFC, P.Value)

x <- subhub2 %>% select(module,moduleMembership) %>% cbind(sub,sub1,sub2)
write_csv(x, "manuscript/brain/manuscript70/results/MEA_black.csv")

