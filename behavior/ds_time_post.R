
# libraries and data
library(tidyverse)

xf <- read_csv("manuscript/behavior/Post_WL.csv")
head(xf)
xf <- xf[-407,]   #batch 10, typo of ids.

# start times
zt <- read_csv("manuscript/behavior/Post_starttimes.csv")
zts <- zt[zt$times=="startTime",]
xf <- full_join(xf, zts %>% select(batch, starttime = ztime))

xf %>% 
  group_by(batch) %>% 
  filter(date == min(date)) %>% 
  filter(ztime < (4200+starttime)) -> xf1

table(xf$group)

# reorganized only
xf1r <- xf1 %>% filter(group=="reorganized")

xf1r

# create list of dataframes for each individual cage
l <- split(xf1r, xf1r$post_batchcage)


# just keep winner loser columns
l.wl <- l %>% map(~ ungroup(.)) %>% map(~ select(., winner,loser) %>% as.data.frame)

#check DS
lapply(l.wl, function(x) compete::ds(compete::get_wl_matrix(x)))

## function to get DS across time for each cage
ds_time <- function(x){
results <- NULL
for(i in 1:nrow(x)){
  results[[i]] <- compete::ds(compete::get_wl_matrix(x[1:i,1:2])) # rows 1 thru i
}
ids <- unlist(lapply(results,names))
vals <- unlist(results)
dx <- data.frame(ids,vals, time = rep(1:nrow(x), times =     unlist(lapply(results, length)) ) )
return(dx)     
}

ds_time(l.wl[[1]])

# just to look
ggplot(ds_time(l.wl[[1]]), aes(x=time, y=vals, color=ids)) + geom_line(lwd=1) + theme_minimal()
ggplot(ds_time(l.wl[[2]]), aes(x=time, y=vals, color=ids)) + geom_line(lwd=1) + theme_minimal()
ggplot(ds_time(l.wl[[3]]), aes(x=time, y=vals, color=ids)) + geom_line(lwd=1) + theme_minimal()


ds_time(l.wl[[23]])

## do for all cages...  
results1 <- NULL
for(i in 1:length(l.wl)){
results1[[i]] <-  ds_time(l.wl[[i]])
}

results1

### add in ranks... 

rank_ds <- function(x){
  x %>% filter(time==max(time))  %>% mutate(rank = rank(-vals, ties.method = 'first')) %>%
    select(ids,rank) %>% full_join(x)
}

results2 <- results1 %>% map(rank_ds)


# but remember for two cages, animals are rank 3/4... (check they are 25hrs too)

# Batch 4 Cage 4  2nd position  [32]  animal id = 2-4
# Batch 2 Cage 3  2nd position  [27]  animal id = 1-4

results2[[27]]$rank <- ifelse(results2[[27]]$rank==1, 1, results2[[27]]$rank +1)
results2[[32]]$rank <- ifelse(results2[[32]]$rank==1, 1, results2[[32]]$rank +1)
results2[[32]]  


## add in batch/cage id, and if alpha/beta/gamma/delta cage...

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

results2 <- Map(cbind, results2, type = c("alphas", "betas", "gammas", "deltas")[as.numeric(substrRight(names(l.wl),1))])
results2 <- Map(cbind, results2, cage=names(l.wl))

results2


## Make one dataframe
DF <- data.table::rbindlist(results2)



## get CIs over time....

DF %>% 
  group_by(rank,time,type) %>%
  summarise(median = median(vals),
            lqr = quantile(vals,.25),
            uqr = quantile(vals,.75)
            ) -> DF.summary


DF.summary

ggplot(DF.summary, aes(x=time,y=median,color=factor(rank))) +
  geom_line() +
  facet_wrap(~type) +
  geom_ribbon(aes(ymin = lqr, ymax = uqr, fill = factor(rank)),alpha=.1, color = NA) +
  theme_classic() +
  ylab("Median David's Scores") +
  xlab("Behavioral Event") +
  scale_fill_manual(values = c("#404788FF", "#238A8DFF", "#55C667FF", "#FDE725FF"),
                      name="Rank", labels=c("1","2","3","4")) +
  scale_color_manual(values = c("#404788FF", "#238A8DFF", "#55C667FF", "#FDE725FF"),
                       name="Rank", labels=c("1","2","3","4")) +
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
        strip.background = element_blank() 
  ) +
  ggtitle("Reemergence of Ranks after Social Reorganization")





### Step graph - latency to win 1,2,3,N fights....

head(xf)

l %>% 
  map(~ ungroup(.)) %>% 
  map(~ select(., post_idbatchcage,ztime,starttime )) %>% 
  map(~ mutate(., ttime = ztime-starttime, contest = row_number())) -> ll


ll <- Map(cbind, ll, type = c("alphas", "betas", "gammas", "deltas")[as.numeric(substrRight(names(l),1))])


ll.df <- data.table::rbindlist(ll)

ggplot(ll.df, aes(x=ttime/60, y=contest, group=post_idbatchcage)) +
  geom_step()+
  facet_wrap(~type)


ll.df %>% group_by(type, contest) %>%
  summarise(median = median(ttime),
            lqr = quantile(ttime,.25),
            uqr = quantile(ttime,.75)
            ) -> ll.df.sum

range(ll.df.sum$contest)

ggplot() +
  geom_step(data=ll.df, aes(x=ttime/60, y=contest, group=post_idbatchcage),alpha=.4)+
  facet_wrap(~type) +
  theme_classic()+
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
       strip.background = element_blank() 
  ) +
  xlab("Latency (minutes") +
  ylab("Contest number") +
  ggtitle("")


## Total aggression boxplots
ll.df %>% separate(post_idbatchcage, sep=" ", into = c("b","batch","c","cageid")) %>%
  mutate(batch_cageid = paste(batch,cageid,sep="-")) %>%
  group_by(batch_cageid,type,batch) %>%
  summarise(total = n()) -> ll.df.total

ll.df$type <- factor(ll.df$type, levels=c("alphas","betas","gammas","deltas"))

library(viridis)

source("functions/geom_boxjitter.R")



ggplot(ll.df.total,aes(x=type,y=total,color=type,fill= type))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4,
                 jitter.height = 0.02, jitter.width = 0.025, errorbar.draw = TRUE,
                 position = position_dodge(0.85))+
  ylab("Total Fights per Cage") +
  xlab("Previous Social Rank") +
  scale_color_manual(values = viridis::viridis(4)) +
  scale_fill_manual(values = viridis::viridis(4))+
  theme_classic() +
  theme(legend.position = 'none')





## Latency boxplots 1st, 5th, 10th...

ll.df %>% separate(post_idbatchcage, sep=" ", into = c("b","batch","c","cageid")) %>%
  mutate(batch_cageid = paste(batch,cageid,sep="-")) %>%
  filter(contest %in% c(1,5,10,20)) ->ll.df.latencies

ll.df.latencies$contest1 <- ifelse(ll.df.latencies$contest==1, "1st", 
                                   ifelse(ll.df.latencies$contest==5, "5th",
                                          ifelse(ll.df.latencies$contest==10, "10th","20th")))

ll.df.latencies$contest1 <- factor(ll.df.latencies$contest1, levels=c("1st","5th","10th","20th"))

ggplot(ll.df.latencies,aes(x=type,y=ttime,color=type,fill= type))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4,
                 jitter.height = 0.02, jitter.width = 0.025, errorbar.draw = TRUE,
                 position = position_dodge(0.85))+
  facet_wrap(~contest1, nrow = 1) +
  ylab("Latency (s)") +
  xlab("Previous Social Rank") +
  scale_color_manual(values = viridis::viridis(4)) +
  scale_fill_manual(values = viridis::viridis(4))+
  theme_classic() +
  theme(axis.text.x = element_text(vjust = 1),
        legend.position = 'none',
        axis.text.y = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.text = element_text(color="#3C3C3C", size=rel(0.8)),
        strip.background = element_blank() 
  )

#mixed models

library(lme4)
library(lmerTest)

ll.df.latencies$type <- factor(ll.df.latencies$type, levels=c("alphas","betas","gammas","deltas"))


ll.df.latencies1 <- ll.df.latencies[ll.df.latencies$contest==1,]
ll.df.latencies5 <- ll.df.latencies[ll.df.latencies$contest==5,]
ll.df.latencies10 <- ll.df.latencies[ll.df.latencies$contest==10,]
ll.df.latencies20 <- ll.df.latencies[ll.df.latencies$contest==20,]
head(ll.df.latencies1)

## 1st 
mixed.lmer1 <- lmer(ttime ~ type + (1|batch), data = ll.df.latencies1)
summary(mixed.lmer1)

plot(mixed.lmer1)  #
qqnorm(resid(mixed.lmer1))
qqline(resid(mixed.lmer1)) 


## 5th 
mixed.lmer5 <- lmer(ttime ~ type + (1|batch), data = ll.df.latencies5)
summary(mixed.lmer5)

plot(mixed.lmer5)  #
qqnorm(resid(mixed.lmer5))
qqline(resid(mixed.lmer5)) 



## 10th
mixed.lmer10 <- lmer(ttime ~ type + (1|batch), data = ll.df.latencies10)
summary(mixed.lmer10)

plot(mixed.lmer10)  #
qqnorm(resid(mixed.lmer10))
qqline(resid(mixed.lmer10)) 


## 20th
mixed.lmer20 <- lmer(ttime ~ type + (1|batch), data = ll.df.latencies20)
summary(mixed.lmer20)

plot(mixed.lmer20)  #
qqnorm(resid(mixed.lmer20))
qqline(resid(mixed.lmer20)) 


## reorder to compare betas to others

ll.df.latencies1$type <- factor(ll.df.latencies1$type, levels=c("betas","gammas","deltas","alphas"))
ll.df.latencies5$type <- factor(ll.df.latencies5$type, levels=c("betas","gammas","deltas","alphas"))
ll.df.latencies10$type <- factor(ll.df.latencies10$type, levels=c("betas","gammas","deltas","alphas"))
ll.df.latencies20$type <- factor(ll.df.latencies20$type, levels=c("betas","gammas","deltas","alphas"))


mixed.lmer1a <- lmer(ttime ~ type + (1|batch), data = ll.df.latencies1)
summary(mixed.lmer1a)

mixed.lmer5a <- lmer(ttime ~ type + (1|batch), data = ll.df.latencies5)
summary(mixed.lmer5a)

mixed.lmer10a <- lmer(ttime ~ type + (1|batch), data = ll.df.latencies10)
summary(mixed.lmer10a)

mixed.lmer20a <- lmer(ttime ~ type + (1|batch), data = ll.df.latencies20)
summary(mixed.lmer20a)


## reorder to compare gammas and deltas

ll.df.latencies1$type <- factor(ll.df.latencies1$type, levels=c("gammas","deltas","alphas","betas"))
ll.df.latencies5$type <- factor(ll.df.latencies5$type, levels=c("gammas","deltas","alphas","betas"))
ll.df.latencies10$type <- factor(ll.df.latencies10$type, levels=c("gammas","deltas","alphas","betas"))
ll.df.latencies20$type <- factor(ll.df.latencies20$type, levels=c("gammas","deltas","alphas","betas"))


mixed.lmer1b <- lmer(ttime ~ type + (1|batch), data = ll.df.latencies1)
summary(mixed.lmer1b)

mixed.lmer5b <- lmer(ttime ~ type + (1|batch), data = ll.df.latencies5)
summary(mixed.lmer5b)

mixed.lmer10c <- lmer(ttime ~ type + (1|batch), data = ll.df.latencies10)
summary(mixed.lmer10c)

mixed.lmer20d <- lmer(ttime ~ type + (1|batch), data = ll.df.latencies20)
summary(mixed.lmer20d)

