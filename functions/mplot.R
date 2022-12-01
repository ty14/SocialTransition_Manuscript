mplot <- function(m, mylevs=NULL, lowcolor="white",highcolor="red1"){
  
  library(ggplot2)
  
  #make the df we will use for plotting
  m.dat <- reshape2::melt(m)
  m.dat <- data.frame(m.dat)
  m.dat <- m.dat[complete.cases(m.dat),] #removing NAs
  
  if(is.null(mylevs)) { mylevs = rownames(m)}
  
  #reorder the levels of the y-axis so plots properly
  m.dat$loser<- factor(m.dat$loser, levels=mylevs)
  m.dat$winner <- factor(m.dat$winner, levels = rev(mylevs))
  m.dat[m.dat == 0] <- NA
  
  #plot
  p1<-ggplot(m.dat, aes(loser, winner, fill = value)) + 
    geom_tile(colour="black", 
              size=0.5, stat="identity") + 
    geom_text(data=m.dat, aes(loser, winner, label = value), color="black", size=rel(3.5)) +
    scale_fill_gradient(low = lowcolor, high = highcolor, space = "Lab", na.value = "white", guide = "colourbar") +
    scale_x_discrete(expand = c(0, 0), position = "top") +
    scale_y_discrete(expand = c(0, 0)) +
    xlab("Loser") + 
    ylab("Winner") +
    theme(axis.text.x = element_text(vjust = 1),
          axis.text.y = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill="white"),
          plot.background = element_rect(fill="white"),
          axis.text = element_text(color="#3C3C3C", size=rel(1.1)),
          legend.position = "none"        
    ) 
  return(p1)
}
