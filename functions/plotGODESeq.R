## Replaces the function from the package SDMTool no longer maintained
##
### ORIGINAL DOCUMENTATION
#' Legend Gradient
#' 
#' \code{legend.gradient} creates and displays a gradient legend on a plot or
#' image file. The place and size of the legend is defined by coordinates,
#' previously identified.
#' 
#' 
#' @param pnts x and y coordinates of the gradient location in the plot
#' @param cols a set of 2 or more colors used in the image, to create the
#' gradient
#' @param limits to label the min and max values of the gradient in the legend
#' @param title to specify the title of the legend
#' @param ... other graphical parameters defined by image() or plot()
#' @return nothing is returned, a gradient legend is added to a plot or a
#' image.
#' @author Lorena Falconi \email{lorefalconi@@gmail.com}
#' @examples
#' 
#' 
#' #define a simple binary matrix
#' tmat = { matrix(c( 0,0,0,1,0,0,1,1,0,1,
#'                    0,0,1,0,1,0,0,0,0,0,
#'                    0,1,NA,1,0,1,0,0,0,1,
#'                    1,0,1,1,1,0,1,0,0,1,
#'                    0,1,0,1,0,1,0,0,0,1,
#'                    0,0,1,0,1,0,0,1,1,0,
#'                    1,0,0,1,0,0,1,0,0,0,
#'                    0,1,0,0,0,1,0,NA,NA,NA,
#'                    0,0,1,1,1,0,0,NA,NA,NA,
#'                    1,1,1,0,0,0,0,NA,NA,NA),nr=10,byrow=TRUE) }
#' 							
#' #do the connected component labeling
#' tasc = ConnCompLabel(tmat)
#' 
#' # Create a color ramp
#' colormap=c("grey","yellow","yellowgreen","olivedrab1","lightblue4")
#'                                                                   
#' #create an image
#' image(tasc,col=colormap, axes=FALSE, xlab="", ylab="", ann=FALSE)
#' 
#' #points for the gradient legend
#' pnts = cbind(x =c(0.8,0.9,0.9,0.8), y =c(1.0,1.0,0.8,0.8))
#' 
#' #create the gradient legend
#' legend.gradient(pnts,colormap,c("Low","High"))
#' 
#' 
#' @export legend.gradient
#' 
##### END ORIGINAL DOCUMENTATION

legend.gradient = function(pnts,cols=heat.colors(100),limits=c(0,1), title='Legend', ...){
  pnts = try(as.matrix(pnts),silent=T)
  if(!is.matrix(pnts)) stop("you must have a 4x2 matrix")
  if(dim(pnts)[1]!=4 || dim (pnts)[2]!=2) stop ("Matrix must have dimensions of 4 rows and 2 columms")
  if(length(cols)<2) stop("You must have 2 or more colors")
  #break up the min and max into a number of values == length(cols)
  yvals = seq(min(pnts[, 2]), max(pnts[, 2]), length=length(cols)+1)
  #cycle through each of the yvals and create polygons
  for (i in 1:length(cols)){  #create the polygon for that color
    polygon(x=pnts[,1],y=c(yvals[i],yvals[i],yvals[i+1],yvals[i+1]),col=cols[i],border=F)
  }
  #add the text
  text(max(pnts[,1]),min(pnts[,2]),labels=limits[1],pos=4,...)
  text(max(pnts[,1]),max(pnts[,2]),labels=limits[2],pos=4,...)
  text(min(pnts[,1]),max(pnts[,2]),labels=title,adj=c(0,-1),...)
}


plotGODESeq <- function(goenrich,
                        deseq, 
                        maxFDR,
                        collapse,
                        color,
                        lowCol,
                        midCol,
                        highCol,
                        scale,
                        label,
                        maxFDRLab,
                        minZscoreLab,
                        extrawidth,
                        centered,
                        fixed_ymax,
                        leghoffset,
                        legvoffset,
                        wrap){
  
  ###### Input:
  
  ###    
  # label: label of the bubbles. Can take values of "id" (default) or "description"
  # char
  # NB: Needs to be checked before goenrich, because value is used in checking goenrich
  if( missing(label) ){ 
    label = "id"
    cat("No label scheme was specified. GO term IDs will be used.\n")
  }
  if( !( label %in% c("id","description") ) ){ 
    label = "id"
    cat("The specified label scheme was neither \"id\" nor \"description\". GO IDs will be used.\n")
  }   
  
  ###  
  # goenrich: result of over-representation analysis. The format is based on 
  # dataframe(
  # no rownames
  # ID          : GO term ID (Factor)
  # description : GO term description (character. Hopefully without ",")
  # Enrich      : ratio between observed and expected annotated genes (numeric)
  # PValue      : significance of the over-representation, from hypergeometric test (numeric) NOT STRICTLY NEEDED AT THE MOMENT
  # FDR         : False Discovery Rate, a.k.a adjusted p-value, normally obtained using Benjamini-Hochberg correction (numeric)
  # genes       : list of genes annotated by the GO term, separated with semicolumns (Factor)
  # )
  # What to do if no ORA input?
  if( missing(goenrich) ){
    stop("No Gene Ontology enrichment dataset provided.")
  }
  if( label == "id" & !( "ID" %in% colnames(goenrich) ) ){
    stop("Gene Ontology enrichment dataset does not provide term IDs, necessary to label the bubbles.")
  }
  if( label == "description" & !( "description" %in% colnames(goenrich) ) ){
    stop("Gene Ontology enrichment dataset does not provide term descriptions, necessary to label the bubbles.")
  }
  if( !( "Enrich" %in% colnames(goenrich) ) ){
    stop("Gene Ontology enrichment dataset does not provide a column \"Enrich\" with term enrichments, necessary for bubble size.")
  }
  # THE FOLLOWING CODE IS NOT NEEDED AT THE MOMENT. KEPT FOR FUTURE USE  
  #  if( !( "PValue" %in% colnames(goenrich) ) ){
  #    stop("Gene Ontology enrichment dataset does not provide non corrected significance of term enrichment, necessary for ???.")
  #  }
  # FDR
  if( !( "FDR" %in% colnames(goenrich) ) ){
    stop("Gene Ontology enrichment dataset does not provide a column \"FDR\" with False Discovery Rates (a.k.a. adjusted P-values). They are necessary for plot Y axis.")
  }
  # genes  
  if( !( "genes" %in% colnames(goenrich) ) ){
    stop("Gene Ontology enrichment dataset does not provide a column \"genes\" with the list of genes annotated by each GO term. This is necessary for computing the zscores used as X axis, and the mean differential expression, used to color bubbles.")
  }
  
  ###    
  # color: color of the bubbles. Can take values of "zscore" (default) or "l2fc"
  # char
  # NB: Needs to be checked before deseq, because its value is used when checking deseq
  if( missing(color) ){ 
    color = "zscore"
    cat("No color scheme was specified. Zscores will be used.\n")
  }
  if( !( color %in% c("zscore","l2fc") ) ){ 
    color = "zscore"
    cat("The specified color scheme was neither \"zscore\" nor \"l2fc\". Zscores will be used.\n")
  }  
  
  ###
  # deseq: result of DESeq2 analysis. We only use two type of information
  # dataframe(
  # rownames: Gene symbols
  # log2FoldChange: Extent of the differential expression for each gene
  # )
  # What to do if no DESeq input?  
  if( missing(deseq) ){
    stop("No differential expression dataset provided.")
  }
  # L2FC 
  if( color == "l2fc" & !( "log2FoldChange" %in% colnames(deseq) ) ){
    stop("Differential Expression dataset does not provide a column \"log2FoldChange\" necessary to compute the plot color.")
  }
  
  ###  
  # maxFDR: maximum value of the ORA False Discovery Rate that we consider for the plot 
  # numeric
  if( missing(maxFDR) ){ 
    maxFDR = 0.05 
    cat(paste("No maximum FDR provided for plotting. Default of", maxFDR,"will be used.\n"))
  }  
  
  ###
  # collapse: either FALSE or value between 0 and 1, representing the proportion of genes in common for two GO terms to be merged. 
  if( missing(collapse) ){ collapse = FALSE  
  cat("No maximum GO term collapsing level is provided. No collapse will be performed.\n")
  } 
  
  ###
  # lowCol: color of the bubbles for the low "color" values
  # valid values for colors, i.e. "blue", "#aabbcc" etc.
  if( missing(lowCol) ){ lowCol = "darkcyan"}
  
  ###
  # midCol: color of the bubbles for the medium "color" values
  # valid values for colors, i.e. "blue", "#aabbcc" etc.
  if( missing(midCol) ){ midCol = "yellow"}
  
  ###
  # highCol: color of the bubbles for the high "color" values
  # valid values for colors, i.e. "blue", "#aabbcc" etc.
  if( missing(highCol) ){ highCol = "darkred"}
  
  ###
  # scale: scale the radius the bubbles.
  # numeric 
  if( missing(scale) ){ scale = 1 }
  
  ###
  # maxFDRLab: maximum Padj value at which bubbles are labelled  
  # numeric
  if( missing(maxFDRLab) ){ 
    maxFDRLab = maxFDR
    cat(paste("No maximum FDR provided for labelling. Default of",maxFDRLab,"will be used.\n"))
  }
  
  ###
  # minZscoreLab: minimum absolute zscore value at which bubbles are labelled  
  # numeric
  if( missing(minZscoreLab) ){ 
    minZscoreLab = 0 
    cat(paste("No minimum zscore provided for labelling. Default of",minZscoreLab,"will be used.\n"))
  }
  
  ###
  # extrawidth: width added to the left of min(zscore) and right of max(zscore)
  # if the "centered" is set to TRUE, it is added on both side to max(abs(min(zscore)),abs(max(zscore)))
  # numeric
  if( missing(extrawidth) ){ extrawidth = 1 }  
  
  ###
  # centered: decides if the plot is symmetrical  
  # Boolean
  if ( missing(centered) ){ centered = FALSE }
  
  ###
  # centered: decides if the plot is symmetrical  
  # Boolean
  if ( missing(fixed_ymax) ){ fixed_ymax = -1 }
  
  ### 
  # value added to the horizontal position of the legend 
  # numeric  
  if( missing(leghoffset)){leghoffset = 0 }
  ###
  # value added to the vertical position of the legend
  # numeric   
  if( missing(legvoffset)){ legvoffset = 0 }
  
  ###
  # wrap: width of label's lines
  # numeric
  if ( missing(wrap) ){ wrap = 15 }  
  
  ##############################
  # loaging required packages
  ##############################
  
  library(GOplot)   # necessary for the function collapsing GO terms that overlap
  library(tidyr)    # necessary for the function separate_rows
  #  library(SDMTools) # necessary for the function legend.gradient
  # Not maintained anymore. The function is now included in this script
  
  ##############################
  ## Preparation of the data
  ##############################
  
  # Replace the p-values and FDR that are zero so that they can be plotted on log-scale.
  
  #  # PVALUE IS NOT USED AT THE MOMENT. KEEP FOR FUTURE USE
  #  minpval <- min(goenrich[goenrich$PValue>0,]$PValue)
  #  # replace the 0 FDR by a tenth of the minimal ones
  #  goenrich$PValue[which(goenrich$PValue==0)] <- minpval/10
  
  if (nrow(goenrich) == 0){
    stop("No GO term in the enrichment file.")
  } else {
    print(paste("Nb of GO terms in input file: ",nrow(goenrich)))
  }
  
  # Get the minimal non-0 FDR
  minfdr <- min(goenrich[goenrich$FDR>0,]$FDR)
  # replace the 0 FDR by a tenth of the minimal ones
  goenrich$FDR[which(goenrich$FDR==0)] <- minfdr/10
  
  # Prepare a subset of the GO enrichment. NB always use this subset, even if it is 100% of the entire results.
  goenrich_subset <- subset(goenrich, goenrich$FDR < maxFDR)
  if (nrow(goenrich_subset) == 0){
    stop("No remaining GO terms in the enrichment file. Try increasing the value of the maxFDR argument.")
  } else {
    print(paste("Nb of GO terms after FDR threshold: ",nrow(goenrich_subset)))
  }
  
  ### Compute zscores, i.e. relative over or underexpression of Genes annotated by each term.
  
  getZscore <- function(goterm) {
    # get the list off all genes annotated by a term
    genes<-strsplit(as.character(unlist(goterm)),";")
    
    # using the DESeq results, get the number of genes that are upregulated or downregulated
    up <- nrow(deseq_data[rownames(deseq_data) %in% unlist(genes)&deseq$log2FoldChange>=0,])
    down <- nrow(deseq_data[rownames(deseq_data) %in% unlist(genes)&deseq$log2FoldChange<0,])
    
    # compute the Zscore
    zscore <- (up - down)/sqrt(up+down)
  }
  
  # Get the Zscore for each GO term
  resZscore<-lapply(goenrich_subset$genes,getZscore)
  
  # add zscore to the GO enrichment table
  goenrich_subset$zscore <- unlist(resZscore)
  
  ### Compute average gene expression enrichment
  
  getMeanL2FC <- function(goterm) {
    # get the list off all genes annotated by a term
    genes<-unlist(strsplit(as.character(unlist(goterm)),";"))
    
    checkGeneDESeq <- function(gene){
      if ( !( gene %in% rownames(deseq_data) ) ){
        cat(paste("gene",gene,"is in the GO enrichment table but not in the list of differential expression\n"))
      }
    }
    lapply(genes,checkGeneDESeq)
    
    # compute the mean of L2FC for all the genes annotated by the term
    meanL2FC <- mean(deseq_data[rownames(deseq_data) %in% genes,]$log2FoldChange)
  }
  
  # Get the mean L2FC for each GO term 
  resL2FC<-lapply(goenrich_subset$genes,getMeanL2FC)
  
  # add average gene expression enrichment to the GO enrichment table
  goenrich_subset$meanL2FC <- unlist(resL2FC)
  
  # rename the FDR into adj_pval and description into term to be compatible with the GOPlot package
  colnames(goenrich_subset)[colnames(goenrich_subset) %in% c("description","FDR")] <- c("term","adj_pval")
  
  # explodes each GO term row into as many rows as genes involved
  # This is necessary for collapsing terms using the package GOPlot
  # separate_rows is a function of the package tidyr
  goenrich_expand <- separate_rows(goenrich_subset, genes)  
  
  ##############################
  ## Plot the data
  ##############################
  
  ## Use GOPlot package to merge together terms annotating similar genesets
  if ( !(collapse == FALSE) ){
    enrich_red <- reduce_overlap(goenrich_expand, overlap = collapse)
    print(paste("Nb of GO terms after collapsing): ",nrow(enrich_red)))
  } else { enrich_red = goenrich_subset } # goenrich_expand is only useful is we collapse.
  
  
  ####### Prepare colours
  
  # build the low palette
  rcPallow <- colorRampPalette(c(lowCol,midCol))
  # build the high palette
  rcPalhigh <- colorRampPalette(c(midCol,highCol))
  
  # FIXME: DEAL WITH THE CASE WHEN ALL THE ZSCORE OR L2FC ARE EITHER POSITIVE OR NEGATIVE
  
  ####### Colors based on zscores
  if (color == "zscore"){
    print(paste("Chosen color for the bubble is: ",color))
    
    # build df to contain all zscores and colors
    datcol <- data.frame(zscore = enrich_red$zscore,colour=NA)
    
    # Extract zscores under 0
    if ( (nrow(enrich_red[enrich_red$zscore < 0,]) == 0) ){
      cat("No negative zscore, only positive color scale is used.\n")
      downreg <- 0
    } else {
      downreg <- enrich_red[enrich_red$zscore < 0,]$zscore
    }
    
    # Extract zscores above 0
    if ( (nrow(enrich_red[enrich_red$zscore > 0,]) == 0) ){
      cat("No strictly positive zscore, only negative color scale is used.\n")
      upreg <- 0
    } else {
      upreg <- enrich_red[enrich_red$zscore >= 0,]$zscore
    }  
    
    # Compute the number of breaks
    maxzscore <- max(abs(downreg),abs(upreg))
    
    if ( !(nrow(enrich_red[enrich_red$zscore < 0,]) == 0) ){
      # I want 10 colors max on each side
      ndownbreak = ceiling(10* max(abs(downreg))/maxzscore)
      # Build negative part of color arrays using the zscore.
      datcollow <- rcPallow(ndownbreak)[as.numeric(cut(downreg,breaks = ndownbreak))]
      # Populate the colours
      datcol[datcol$zscore < 0,]$colour <- datcollow
    }
    
    if ( !(nrow(enrich_red[enrich_red$zscore > 0,]) == 0) ){
      # I want 10 colors max on each side
      nupbreak = ceiling(10* max(abs(upreg))/maxzscore  )
      # Build positive part of color arrays using the zscore.
      datcolhigh <- rcPalhigh(nupbreak)[as.numeric(cut(upreg,breaks = nupbreak))]
      # Populate the colours
      datcol[datcol$zscore >= 0,]$colour <- datcolhigh
    }
  }
  
  ####### Colors based on L2FC
  if (color == "l2fc"){
    print(paste("Chosen color for the bubble is: ",color))
    
    # build df to contain all L2FC and colors
    datcol <- data.frame(L2FC = enrich_red$meanL2FC,colour=NA)
    
    # Extract L2FC under 0
    if ( (nrow(enrich_red[enrich_red$meanL2FC < 0,]) == 0) ){
      cat("No negative meanL2FC, only positive color scale is used.\n")
      downreg <- 0
    } else {
      downreg <- enrich_red[enrich_red$meanL2FC < 0,]$meanL2FC
    }
    
    # Extract L2FC above 0
    if ( (nrow(enrich_red[enrich_red$meanL2FC > 0,]) == 0) ){
      cat("No strictly positive L2FC, only negative color scale is used.\n")
      upreg <- 0
    } else {
      upreg <- enrich_red[enrich_red$meanL2FC >= 0,]$meanL2FC
    }
    
    # Compute the number of breaks
    maxl2fc <- max(abs(downreg),abs(upreg))
    
    if ( !(nrow(enrich_red[enrich_red$meanL2FC < 0,]) == 0) ){
      # I want 10 colors max on each side
      ndownbreak = ceiling(10* max(abs(downreg))/maxl2fc)
      # Build negative part of color arrays using the L2FC
      datcollow <- rcPallow(ndownbreak)[as.numeric(cut(downreg,breaks = ndownbreak))]
      # Populate the colours
      datcol[datcol$L2FC < 0,]$colour <- datcollow
    }
    
    if ( !(nrow(enrich_red[enrich_red$meanL2FC > 0,]) == 0) ){
      # I want 10 colors max on each side
      nupbreak = ceiling(10* max(abs(upreg))/maxl2fc  )
      # Build positive part of color arrays using the zscore.
      datcolhigh <- rcPalhigh(nupbreak)[as.numeric(cut(upreg,breaks = nupbreak))]
      # Populate the colours
      datcol[datcol$L2FC >= 0,]$colour <- datcolhigh
    }
    
  }
  
  ####### create width and height of plot
  
  if(centered == TRUE){
    xmin = -max(abs(min(enrich_red$zscore)),abs(max(enrich_red$zscore))) - extrawidth
    xmax = max(abs(min(enrich_red$zscore)),abs(max(enrich_red$zscore))) + extrawidth
  } else {
    xmin = min(enrich_red$zscore) - extrawidth
    xmax = max(enrich_red$zscore) + extrawidth
  }
  ymin = min(-log10(enrich_red$adj_pval)) - 0.5
  if(fixed_ymax != -1){
    cat("Manual ymax entered: ",fixed_ymax,"\n")
    ymax = fixed_ymax
  } else {
    ymax = max(-log10(enrich_red$adj_pval)) + 0.5
  }
  
  ###### Plot!
  
  # FIXME: PROVIDE AN OPTION TO USE THE RADIUS INSTEAD OF SURFACE FOR THE ENRICHMENT
  symbols(enrich_red$zscore, 
          -log10(enrich_red$adj_pval),
          circle = sqrt(enrich_red$Enrich/pi)*scale, # gives the surface in fonction of enrichment. 
          inches=FALSE,
          fg="white",
          bg = datcol$colour,
          xlim=c(xmin,xmax),ylim=c(ymin,ymax),
          xlab="zscore",ylab="-log(adjusted p-value GO enrichment)"
  )
  ###### Create the Legend 
  
  # Create a color ramp
  rbPal <- colorRampPalette(c(lowCol,midCol,highCol))
  legcol <- rbPal(50)
  
  #points for the gradient legend
  
  pnts = cbind(x =c(xmin+0.5+leghoffset,xmin+0.75+leghoffset,xmin+0.75+leghoffset,xmin+0.5+leghoffset), 
               y =c(ymin+2.5+legvoffset,ymin+2.5+legvoffset,ymin+1+legvoffset,ymin+1+legvoffset))
  
  #create the gradient legend
  
  if (color == "l2fc") {
    legend.gradient(pnts,
                    legcol,
                    c(sprintf("%.2f",min(enrich_red$meanL2FC)),sprintf("%.2f",max(enrich_red$meanL2FC))), 
                    title = "mean Expr L2FC")
  } else {
    legend.gradient(pnts,
                    legcol,
                    c(sprintf("%.2f",min(enrich_red$zscore)),sprintf("%.2f",max(enrich_red$zscore))), 
                    title = "Zscore")
  }
  
  # Legend for enrichment: black circle of surface "1". 
  # NB: I have to feed a vector of 1 element for the size to work. Don't ask. Symbols is crazy
  symbols(c(xmin+0.5+leghoffset+0.125), # move by 0.125 because gradient is 0.25 large
          c(ymin+3.5+legvoffset),
          circle = c(sqrt(1/pi))*scale, 
          inches=FALSE,
          fg="white",
          bg = "black",
          add = TRUE
  )
  text(xmin+0.5+leghoffset+0.125+sqrt(1/pi)*scale,ymin+3.5+legvoffset,label="Enrich = 1",pos=4)
  
  # Select the bubbles that will be labeled. 
  # We only select adj p-value less than maxFDRLab and zscore inferieur to -minZscoreLab or superieur to minZscoreLab
  tolabel <- subset(enrich_red, enrich_red$adj_pval < maxFDRLab | abs(enrich_red$zscore)>minZscoreLab)
  
  
  if (label == "id"){
    ####### GO ID as labels
    
    text(enrich_red$zscore, 
         -log10(enrich_red$adj_pval),
         enrich_red$ID,
         cex=0.75
    )
  } 
  
  if (label == "description"){
    ####### GO description as labels
    
    ## Build wrapped label
    
    # Core wrapping function
    wrap.it <- function(x, len)
    { 
      sapply(x, function(y) paste(strwrap(y, len), 
                                  collapse = "\n"), 
             USE.NAMES = FALSE)
    }
    
    # Call this function with a list or vector
    wrap.labels <- function(x, len)
    {
      if (is.list(x))
      {
        lapply(x, wrap.it, len)
      } else {
        wrap.it(x, len)
      }
    }
    
    labels <- wrap.labels(tolabel$term, wrap)
    
    text(tolabel$zscore, 
         -log10(tolabel$adj_pval),
         labels,
         cex=0.75
    )
  }
  
} ### END OF MAIN FUNCTION






goenrich_data <- read.table("GO-example2.csv",sep="\t",fill=T,quote="\"",header=T)

# rename the columns to make them less weird and compatible with the GOPlot package.
# ID is the GO ID, like GO:0001234
# Enrich is the enrichment ratio 
# genes is the list of genes contributing to this term enrichment, separated by semi-columns
colnames(goenrich_data)[colnames(goenrich_data) %in% c("geneset","R","OverlapGene_UserID")] <- c("ID","Enrich","genes")

# remove commas from descriptions, because they suck
goenrich_data$description <- gsub(',',"",goenrich_data$description)

# Load results from DESeq2
ViSEAGO::show_table(asc_Results) 
go <- asc_Results@data %>% as.data.frame() %>% select( GO.ID , aes_upx.pvalue)
head(go)
# remove commas from GO term descriptions, because they suck
goenrich_data$description <- gsub(',',"",goenrich_data$description)

deseq_data <- sa_top %>% rename(log2FoldChange = logFC)

plotGODESeq(goenrich_data,
            deseq_data,
            maxFDR = 0.05,
            collapse = 1,
            color="l2fc",
            lowCol = "deepskyblue4",
            midCol = "#DDDDDD",
            highCol = "firebrick",
            extrawidth=3,
            centered=FALSE,
            fixed_ymax = 10,
            leghoffset=-0.5,
            legvoffset=1.5,
            label = "description",
            scale = 0.7,
            maxFDRLab = 1e-12,
            minZscoreLab = 2.5,
            wrap = 15
)
