###  Post-reorganization Behavior Analysis - Making Example Matrix

## This uses batch 7 as an example, but could be extended to any batch.



## Source pre_org_behavior ----

source("manuscript/behavior/pre_org_behavior.R")

library(tidyverse)

xf <- read_csv("manuscript/behavior/Post_WL.csv")
head(xf)



## Some Functions ----

# function to melt matrix but keep ids
melt.mat1 <- function(q){
  q1<-reshape2::melt(q)
  q1$winner1 <- 1:4
  q1$loser1 <- rep(1:4, each =4)
  return(q1)
}


substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


addid <- function(x){
  x$winnerid <- paste(substrRight(as.character(x$names),1),x$winner,sep="-")
  x$loserid <- paste(substrRight(as.character(x$names),1),x$loser,sep="-")
  return(x)
}

# create new ids based on A,B,C,D
addid2 <- function(x){
  x$winner2 <- paste0(x$cageid, x$winner1)
  x$loser2 <- paste0(x$cageid, x$loser1)
  return(x)
}

add_value1 <- function(x) { x$value1 <- ifelse(x$value==0, NA, x$value); return(x)}


#rescale value between two numbers e.g. 0-70 and assign color
add_colors <- function(x, minx = 0, maxx = 70){
  d <- scales::rescale(x %>% map(~.$value) %>% unlist %>% unname, to=c(minx,maxx))
  rr <- range(d)
  svals <- (d-rr[1])/diff(rr)
  f <- colorRamp(c("white", "red1"))
  colors <- rgb(f(svals)/255)
  x <-  Map(cbind, x, value2 = split(d, ceiling(seq_along(d)/16)))
  x <-  Map(cbind, x, clr = split(colors, ceiling(seq_along(colors)/16)))
  return(x)
}

# Matrix plot function
plot_mat0 <- function(x){
  ggplot(x, aes(loser1, winner1, label=value1, fill=clr)) + 
    geom_tile(colour="black", 
              size=0.5, stat="identity") + 
    scale_y_continuous(trans = "reverse", breaks = 1:4, labels = unique(x$winner2)) +
    scale_x_continuous(breaks=1:4, labels = unique(x$loser2)) +
    geom_text(color="black", size=rel(2.5)) +
    scale_fill_identity() +
    #scale_fill_gradient(low = 'white', high = 'red1', space = "Lab", na.value = "white", guide = "colourbar") +
    theme(axis.text.x = element_text(vjust = 1),
          axis.text.y = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          axis.text = element_text(color="#3C3C3C", size=rel(0.8)),
          legend.position = "none",
          strip.background = element_blank(),
          strip.text.x = element_blank()
    )
}


add_ids <- function(x,ids){
  x$winnerid <- gsub("Cage ", "", x$winner)
  x$loserid <- gsub("Cage ", "", x$loser)
  x$winner2 <- ids$winner2[match(x$winnerid, ids$winnerid)]
  x$loser2 <- ids$winner2[match(x$loserid, ids$winnerid)]
  return(x)
}



## Plot one pre-org batch ----


# melt pre matrices - 
# need to order these in the list based on most to least aggression  
# get order of highest to lowest
premat7 <- l.mat[grepl("Batch 7", names(l.mat))] %>% map(melt.mat1) 
premat7 <- Map(cbind, premat7, names = names(premat7))
premat7 <- premat7 %>% map(addid)%>% map(~ mutate(., total = sum(value)))
matorder <-   premat7 %>% map(~ filter(.,row_number()==1)) %>% map(~.$total) %>% unlist() %>% order() %>% rev()
premat7x <- Map(cbind, premat7[matorder], cageid = LETTERS[1:4])
df7pre <- premat7x %>% map(addid2) %>% map(add_value1)
df7pre <- add_colors(df7pre, minx=0, maxx=70) #add colors


# Plots - pre

p1=plot_mat0(df7pre[[1]]) + ggtitle("Cage A") + xlab("") + ylab("")
p2=plot_mat0(df7pre[[2]]) + ggtitle("Cage B") + xlab("") + ylab("")
p3=plot_mat0(df7pre[[3]]) + ggtitle("Cage C") + xlab("") + ylab("")
p4=plot_mat0(df7pre[[4]]) + ggtitle("Cage D") + xlab("") + ylab("")

library(gridExtra)

grid.arrange(p1,p2,p3,p4, nrow=1)




#### Post matrices....

# Melt matrices
xf7 <- xf %>% filter(batch==7) %>% select(Cage, winner, loser)
l7 <- split(xf7, xf7$Cage) # create list of dataframes for each individual cage
l7.wl <- l7 %>% map(~ select(., winner,loser) %>% as.data.frame)  # just keep winner loser columns

# convert list to win-loss matrices organized by David's Scores
l7.mat <- l7.wl %>% map(~ compete::get_wl_matrix(.)) %>% map(~ compete::org_matrix(., method="ds"))

# melt post matrices
ranknames <- c("alphas","betas", "gammas", "deltas")
df7 <- Map(cbind, l7.mat  %>% map(melt.mat1), type = ranknames)


# make pre melted matrices into df
df7prex <- premat7x %>% map(addid2) %>% data.table::rbindlist()
ids <- df7prex %>% select(winnerid, winner2) %>% unique()

df7 <- df7 %>% map(add_ids,ids) %>% map(add_value1)  # add ids
df7 <- add_colors(df7)  # add colors


### Plot matrices

p5=plot_mat0(df7[[1]]) + ggtitle("Alphas") + xlab("") + ylab("")
p6=plot_mat0(df7[[2]]) + ggtitle("Betas") + xlab("") + ylab("")
p7=plot_mat0(df7[[3]]) + ggtitle("Gammas") + xlab("") + ylab("")
p8=plot_mat0(df7[[4]]) + ggtitle("Deltas") + xlab("") + ylab("")

grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,nrow=2)
