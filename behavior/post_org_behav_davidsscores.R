###  Post-reorganization Behavior Analysis



## David's Scores and Directional Consistency within 1 hour of reorganization.


# libraries and data
library(tidyverse)

xf <- read_csv("manuscript/behavior/Post_WL.csv")
head(xf)
xf <- xf[-407,]   #batch 10, typo of ids.




## Just assuming that lowest ztime on first day is start of observations:

xf %>% filter(time=="1 hour") %>% group_by(batch) %>% summarise(diff = max(ztime)-min(ztime))

# get data within 70 minutes (4200 seconds) of first day of reorganization.

xf %>% 
  group_by(batch) %>% 
  filter(date == min(date)) %>% 
  filter(ztime < (4200+min(ztime))) -> xf1

table(xf$group)

# reorganized only
xf1r <- xf1 %>% filter(group=="reorganized")

xf1r

# create list of dataframes for each individual cage
l <- split(xf1r, xf1r$post_batchcage)

# just keep winner loser columns
l.wl <- l %>% map(~ ungroup(.)) %>% map(~ select(., winner,loser) %>% as.data.frame)

# convert list to win-loss matrices organized by David's Scores
l.mat <- l.wl %>% map(~ compete::get_wl_matrix(.)) %>% map(~ compete::org_matrix(., method="ds"))


# Get average David's Scores for each rank.
l.ds <- l.mat %>% map(~ compete::ds(.)) 


## we have to reinsert "0" DS for two animals that don't appear in any wins/losses
#code below helps id and fix this:
l.ds  
which(names(l.ds)=="Batch 4 Cage 4") #[32]
which(names(l.ds)=="Batch 2 Cage 3") #[27]

# Batch 4 Cage 4  2nd position  [32]  animal id = 2-4
# Batch 2 Cage 3  2nd position  [27]  animal id = 1-4

# manually inputting
x = l.ds[32]
x1 <- c(x[[1]][1],0,x[[1]][2:3])
names(x1)[2]<-"Cage 2-4"
l.ds[[32]] <- x1

x = l.ds[27]
x1 <- c(x[[1]][1],0,x[[1]][2:3])
names(x1)[2]<-"Cage 1-4"
l.ds[[27]] <- x1


###
ds.df <- l.ds  %>% map(~ matrix(., ncol=4)) %>% do.call("rbind",.)
colnames(ds.df) <-c("Rank1", "Rank2", "Rank3", "Rank4")
ds.df

## Make boxplot of DS

# Reorganize Data

ds.long <- ds.df %>% as.data.frame %>% pivot_longer(1:4, names_to = "Rank")
ds.long$id <- rep(1:(nrow(ds.long)/4), each=4)

#get summary stats
ds.long %>% 
  group_by(Rank) %>% 
  summarise(median = median(value),
            lqr = quantile(value,.25),
            uqr = quantile(value,.75)
  ) -> ds.long.summary



ggplot() + 
  geom_line(data=ds.long, aes(x=Rank, y=value, group=id), alpha=.3, color="gray57") +
  theme_classic() +
  xlab("Mouse Rank") +
  ylab("David's Scores") +
  geom_hline(yintercept=0, lty=2, color="red", alpha=.5) +
  geom_errorbar(data=ds.long.summary, 
                aes(x=Rank, ymin=lqr, ymax=uqr), width=0.0, size=1, color="firebrick") +
  geom_point(data=ds.long.summary, 
             aes(x=Rank, y=median), size=3, shape=21, fill="white") +
  ggtitle("David's Scores after 70 minutes of Group Reorganization") +
  ylim(-6,6)


### Breakdown by Cage Rank

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

ds.long$cage <- rep(substrRight(names(l.ds),1), each=4)
ds.long$cagerank <- rep(rep(c("Alphas", "Betas", "Gammas", "Deltas"), each = 4), nrow(ds.long)/16)
ds.long$cagerank <- factor(ds.long$cagerank, levels =c("Alphas", "Betas", "Gammas", "Deltas"))


#get summary stats
ds.long %>% 
  group_by(cagerank,Rank) %>% 
  summarise(median = median(value),
            lqr = quantile(value,.25),
            uqr = quantile(value,.75)
  ) -> ds.long.cage.summary


ggplot() + 
  geom_line(data=ds.long, aes(x=Rank, y=value, group=id), alpha=.3, color="black") +
  theme_classic() +
  xlab("New Rank") +
  ylab("David's Scores") +
  geom_hline(yintercept=0, lty=2, color="red", alpha=.5) +
   geom_errorbar(data=ds.long.cage.summary, 
                aes(x=Rank, ymin=lqr, ymax=uqr), width=0.0, size=1, color="firebrick") +
  geom_point(data=ds.long.cage.summary, 
             aes(x=Rank, y=median), size=3, shape=21, fill="white") +
  facet_wrap(~cagerank)+
  ylim(-6,6) +
  theme(axis.text.x = element_text(vjust = 1),
        axis.text.y = element_text(hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 15),
        axis.text = element_text(color="#3C3C3C", size=rel(1)),
        strip.text = element_text(size=11),
        strip.background = element_blank()) 

# +
  ggtitle("David's Scores 70 minutes after Group Reorganization") 

+
### Positive David's Scores

ds.long %>% filter(Rank == "Rank2", value>0)
ds.long %>% filter(Rank == "Rank2", value<=0)

ds.long %>% filter(Rank == "Rank3", value>0)
ds.long %>% filter(Rank == "Rank3", value<=0)


ds.long %>% group_by(Rank,cagerank) %>% summarise(total = sum(value>0)) %>%
  pivot_wider(values_from = total, names_from = cagerank)



## Directional Consistency.

# get results of DC test across all cages
l.dc <- l.mat %>% map(~ compete::dc_test(.))


l.dc.dc <- l.dc %>% map(~ .[[7]])    #DC vals

median(unlist(l.dc.dc))
quantile(unlist(l.dc.dc),.25)
quantile(unlist(l.dc.dc),.75)

l.dc.pval <- l.dc %>% map(~ .[[1]])  # pvalues
l.dc.pval %>% unlist() %>% round(3)


## by ranks

l.mat.alphas <- l.mat[seq(1,52,4)]
l.mat.betas <- l.mat[seq(2,52,4)]
l.mat.gammas <- l.mat[seq(3,52,4)]
l.mat.deltas <- l.mat[seq(4,52,4)]

l.dc.alphas <- l.mat.alphas %>% map(~ compete::dc_test(.))
l.dc.betas <- l.mat.betas %>% map(~ compete::dc_test(.))
l.dc.gammas <- l.mat.gammas %>% map(~ compete::dc_test(.))
l.dc.deltas <- l.mat.deltas %>% map(~ compete::dc_test(.))

l.dc.dc.alphas <- l.dc.alphas %>% map(~ .[[7]])    #DC vals
l.dc.dc.betas <- l.dc.betas %>% map(~ .[[7]])    #DC vals
l.dc.dc.gammas <- l.dc.gammas %>% map(~ .[[7]])    #DC vals
l.dc.dc.deltas <- l.dc.deltas %>% map(~ .[[7]])    #DC vals


median(unlist(l.dc.dc.alphas))
quantile(unlist(l.dc.dc.alphas),.25)
quantile(unlist(l.dc.dc.alphas),.75)

median(unlist(l.dc.dc.betas))
quantile(unlist(l.dc.dc.betas),.25)
quantile(unlist(l.dc.dc.betas),.75)

median(unlist(l.dc.dc.gammas))
quantile(unlist(l.dc.dc.gammas),.25)
quantile(unlist(l.dc.dc.gammas),.75)

median(unlist(l.dc.dc.deltas))
quantile(unlist(l.dc.dc.deltas),.25)
quantile(unlist(l.dc.dc.deltas),.75)

