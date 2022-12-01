
#### Post-Organization Matrices, David's Scores, Transitivity etc.

# load packages and data
library(tidyverse)
library(compete)

df <- read_csv("manuscript/behavior/Post_WL.csv")

head(df)


fix10 <- df %>% filter(post_batchcage == "Batch 10 Cage 3")
fix10$loser
fix10["loser"][fix10["loser"]=="Cage 3-3"]<- "Cage 3-1"

df <- df %>% filter(post_batchcage != "Batch 10 Cage 3")
df <- df %>% rbind(fix10)

## 1. Create aggregated WL matrices for each cage.
# these need to be ordered by David's Scores

# Just examining the number of rows of data for each cage/batch
table(df$batch, df$Cage)
table(df$post_batchcage)



# $`Batch 10 Cage 3` <- 5 mice id Cage 3-3 needs to be changed to Cage 3-1
# create list of dataframes for each individual cage
l <- split(df, df$post_batchcage)



# just keep winner loser columns
l.wl <- l %>% map(~ select(., winner,loser) %>% as.data.frame) 

# convert list to win-loss matrices organized by David's Scores
l.mat <- l.wl %>% map(~ compete::get_wl_matrix(.)) %>% map(~ compete::org_matrix(., method="ds"))


# just work out which are control and which are reorganized
df.groups <- df %>% select(post_idbatchcage, group, time) %>% unique

df.groups


# split l.mat into 'control' and 'reorganized'
ids <- df.groups$group[match(names(l.mat), df.groups$post_idbatchcage)]

l.mat.control <- l.mat[which(ids=="control")]
l.mat.reorg <- l.mat[which(ids=="reorganized")]


# Get average David's Scores for each rank. - for reorg vs control.

ds.reorg <- l.mat.reorg %>% map(~ compete::ds(.)) %>% map(~ matrix(.,ncol = 4)) %>% do.call("rbind",.)

ds.control <- l.mat.control %>% map(~ compete::ds(.)) %>% map(~ matrix(., ncol=4)) %>% do.call("rbind",.)

colnames(ds.control)<-colnames(ds.reorg)<-c("Rank1", "Rank2", "Rank3", "Rank4")



## Make boxplot of DS

# Reorganize Data
ds.reorg <- as.data.frame(ds.reorg)
ds.control <- as.data.frame(ds.control)

ds.reorg$group <- "Reorganized" 
ds.control$group <- "Control" 
ds.reorg$time <- rep(c("25hr", "70min", "25hr","70min"), times= c(4,20,24,4))


ds.reorg.long <- ds.reorg %>% pivot_longer(1:4, names_to = "Rank")
ds.reorg.long$id <- rep(1:(nrow(ds.reorg.long)/4), each=4)

ds.control.long <- ds.control %>% pivot_longer(1:4, names_to = "Rank")
ds.control.long$id <- rep(1:(nrow(ds.control.long)/4), each=4)



#get summary stats
ds.reorg.long %>% 
  group_by(Rank,time) %>% 
  summarise(median = median(value),
            lqr = quantile(value,.25),
            uqr = quantile(value,.75)
  ) -> ds.reorg.summary


ds.control.long %>% 
  group_by(Rank) %>% 
  summarise(median = median(value),
            lqr = quantile(value,.25),
            uqr = quantile(value,.75)
  ) -> ds.control.summary


#combine
ds.all.long <- rbind(ds.reorg.long,ds.control.long) %>% mutate(id = paste0(group,id))

ds.all.long %>% 
  group_by(Rank,time) %>% 
  summarise(median = median(value),
            lqr = quantile(value,.25),
            uqr = quantile(value,.75)
  ) -> ds.all.summary

ds.all.summary

ds.all.long %>% group_by(Rank) %>% count()

ds.all.long %>% filter(Rank == "Rank2", value>0)
ds.all.long %>% filter(Rank == "Rank2", value<=0)

ds.all.long %>% filter(Rank == "Rank3", value>0)
ds.all.long %>% filter(Rank == "Rank3", value<=0)

ds.reorg.long$time <- factor(ds.reorg.long$time, levels = c("70min", "25hr"))
ds.reorg.summary$time <- factor(ds.reorg.summary$time, levels = c("70min", "25hr"))
ds.reorg.long70<- ds.reorg.long %>% filter(time == "70min")
ds.reorg.long25 <-ds.reorg.long %>% filter(time == "25hr")
ds.reorg.summary70<- ds.reorg.summary%>% filter(time == "70min")
ds.reorg.summary25<- ds.reorg.summary%>% filter(time == "25hr")


p1 <- ggplot() + 
  geom_line(data=ds.reorg.long70, aes(x=Rank, y=value, group=id), alpha=.3, color="gray57") +
  theme_classic() +
  xlab("") +
  ylab("David's Scores") +
  geom_hline(yintercept=0, lty=2, color="red", alpha=.5) +
  geom_errorbar(data=ds.reorg.summary70, 
                aes(x=Rank, ymin=lqr, ymax=uqr), width=0.0, size=1, color="firebrick") +
  geom_point(data=ds.reorg.summary70, 
             aes(x=Rank, y=median), size=3, shape=21, fill="white") +
  ggtitle("Reorganized 70 min") +
  ylim(-6,6)+
  theme(axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),    
        axis.title.x = element_text(hjust = 1),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        text = element_text(size = 20)
  
  )



p2 <- ggplot() + 
  geom_line(data=ds.reorg.long25, aes(x=Rank, y=value, group=id), alpha=.3, color="gray57") +
  theme_classic() +
  xlab("") +
  ylab("") +
  geom_hline(yintercept=0, lty=2, color="red", alpha=.5) +
  geom_errorbar(data=ds.reorg.summary25, 
                aes(x=Rank, ymin=lqr, ymax=uqr), width=0.0, size=1, color="firebrick") +
  geom_point(data=ds.reorg.summary25, 
             aes(x=Rank, y=median), size=3, shape=21, fill="white") +
  ggtitle("Reorganized 25 hr") +
  ylim(-6,6)+
  theme(axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),    
        axis.title.x = element_text(hjust = 1),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        text = element_text(size = 20)
        
  )




ds.control.long$time <- factor(ds.control.long$time, levels = c("70min", "25hr"))
ds.control.summary$time <- factor(ds.control.summary$time, levels = c("70min", "25hr"))
p3 <- ggplot() + 
  geom_line(data=ds.control.long, aes(x=Rank, y=value, group=id), alpha=.3, color="gray57") +
  theme_classic() +
  xlab("") +
  ylab("") +
  geom_hline(yintercept=0, lty=2, color="red", alpha=.5) +
  geom_errorbar(data=ds.control.summary, 
                aes(x=Rank, ymin=lqr, ymax=uqr), width=0.0, size=1, color="firebrick") +
  geom_point(data=ds.control.summary, 
             aes(x=Rank, y=median), size=3, shape=21, fill="white") +
  ggtitle("Control") +
  ylim(-6,6)+
  theme(axis.text.x = element_text(vjust = 1, size = 20),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),    
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        text = element_text(size = 20)
        
  )

  

bottom <- grid::textGrob("Mouse Rank", gp = grid::gpar(fontsize = 20))
library(gridExtra)
ds <- grid.arrange(p1,p2,p3, nrow=1,bottom = bottom)

ggsave("manuscript/behavior/results_figures/post_ds2.png",ds,height =6, width =15, dpi=600)


## Directional Consistency

# get results of DC test across all cages
l.dc <- l.mat %>% map(~ compete::dc_test(.))


l.dc.dc <- l.dc %>% map(~ .[[7]])    #DC vals

median(unlist(l.dc.dc))
quantile(unlist(l.dc.dc),.25)
quantile(unlist(l.dc.dc),.75)

l.dc.pval <- l.dc %>% map(~ .[[1]])  # pvalues
l.dc.pval %>% unlist() %>% round(3)




### Image of matrices - broken down by who reorganized with who.

# melt the list of matrices into a dataframe

l.mat.reorg
names(l.mat.reorg)

# get batch names
batches <- trimws(substr(names(l.mat.reorg),1,8), "r")

# each row being a batch.  
# order the matrices in each row by the degree of aggression
# order the rows by degree of aggression

melt_mat <- function(x){
  colnames(x)<-rownames(x)<-1:4
  return(reshape2::melt(x))
}


melt.mat.reorg <- Map(cbind, l.mat.reorg %>% map(melt_mat), Batch = batches, Cageid = 1:length(l.mat.reorg))
melt.mat.reorg <- Map(cbind, melt.mat.reorg, batchcage = names(melt.mat.reorg)) # add in batchcage
melt.mat.reorg.df <- data.table::rbindlist(melt.mat.reorg)

# adding in total fights per cage
melt.mat.reorg.df <- melt.mat.reorg.df %>% group_by(Batch,Cageid) %>% mutate(total = sum(value))

# adding in total fights per batch
melt.mat.reorg.df <- melt.mat.reorg.df %>% group_by(Batch) %>% mutate(totalB = sum(value))

# order of fights by batch highest to lowest
batch.levs <- melt.mat.reorg.df %>% group_by(Batch) %>% filter(row_number()==1) %>% arrange(-totalB) %>% .$Batch

# add in number of aggr level for batches
melt.mat.reorg.df$Batchno <- match(melt.mat.reorg.df$Batch,batch.levs)

# work out which cages in each batch have highest aggression
cageranks <- melt.mat.reorg.df %>% group_by(Batch, Cageid) %>% filter(row_number()==1) %>%
  arrange(-total, .by_group = TRUE) %>% 
  group_by(Batch) %>%
  mutate(cagerank = rank(-total)) %>%
  select(Cageid, cagerank) 

table(cageranks$Batch, cageranks$cagerank) #checking we have one rank per batch

# adding in rank of each cage back into df
melt.mat.reorg.df$cagerank <- cageranks$cagerank[match(melt.mat.reorg.df$Cageid,cageranks$Cageid)]

melt.mat.reorg.df

# make text column with no 0s
melt.mat.reorg.df$value1<-ifelse(melt.mat.reorg.df$value!=0, melt.mat.reorg.df$value, NA)
melt.mat.reorg.df

pm1 <- ggplot(melt.mat.reorg.df, aes(loser, winner, fill = value)) + 
  geom_tile(colour="black", 
            size=0.5, stat="identity") + 
  scale_y_continuous(trans = "reverse") +
  geom_text(data=melt.mat.reorg.df, aes(loser, winner, label = value1), color="black", size=5) +
  scale_fill_gradient(low = 'white', high = 'red1', space = "Lab", na.value = "white", guide = "colourbar") +
  facet_grid(vars(Batchno), vars(cagerank)) +
  theme(axis.text.x = element_text(vjust = 1, size = 20),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),    
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 15)
  )


pm1


## Make Dichotomized Plot

# make binary matrices
l.dimat.reorg <- l.mat.reorg %>% map(~ compete::get_di_matrix(.))

melt.dimat.reorg.df <-
  Map(cbind, l.dimat.reorg %>% map(melt_mat), batchcage = names(l.dimat.reorg)) %>% data.table::rbindlist()

colnames(melt.dimat.reorg.df)[3]<-"divalue"

# join in dichotomized value
melt.mat.reorg.df1 <- full_join(melt.mat.reorg.df, melt.dimat.reorg.df)


# get DC for each relationship
melt.dcmat.reorg.df <-
  Map(cbind,
      l.mat.reorg %>% 
        map(function(x) x / (x + t(x) ) ) %>% 
        map(melt_mat),
      batchcage = names(l.dimat.reorg)) %>% 
  data.table::rbindlist()

colnames(melt.dcmat.reorg.df)[3]<-"dcvalue"

# join in dc value
melt.mat.reorg.df1 <- full_join(melt.mat.reorg.df1, melt.dcmat.reorg.df)

# ensure any 0/0 relationships are NA
melt.mat.reorg.df1$divalue <- ifelse(is.na(melt.mat.reorg.df1$dcvalue), NA, melt.mat.reorg.df1$divalue)

# make all 0s NAs
melt.mat.reorg.df1$divalue <- ifelse(melt.mat.reorg.df1$divalue==0, NA, melt.mat.reorg.df1$divalue)


pm2 <- ggplot(melt.mat.reorg.df1, aes(loser, winner, fill = dcvalue)) + 
  geom_tile(colour="black", 
            size=0.5, stat="identity") + 
  scale_y_continuous(trans = "reverse") +
  geom_text(data=melt.mat.reorg.df1, aes(loser, winner, label = divalue), color="black", size=5) +
  scale_fill_gradient(low = 'white', high = 'red1', space = "Lab", na.value = "white", guide = "colourbar") +
  facet_grid(vars(Batchno), vars(cagerank)) +
  theme(axis.text.x = element_text(vjust = 1, size = 20),
        axis.text.y = element_text(hjust = 0.5, size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.text = element_text(color="#3C3C3C", size = 20),
        axis.title = element_text(size = 20),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 15)
        
  )


pm2

pm3 <- grid.arrange(pm1,pm2,nrow=1)
ggsave("manuscript/behavior/results_figures/post_matrics.png",pm3,height =16, width =14, dpi=600)

# checking which batch was row 1 of matrix
melt.mat.reorg.df %>% filter(Batchno==1)

l.mat[grepl("Batch 7", names(l.mat))]
