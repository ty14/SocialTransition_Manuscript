## cort and aggression analysis 

#libraries and data 
library(tidyverse)
library(compete)

# functions 
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


# First lets get aggression per cage post 
ax <- read_csv("manuscript/behavior/Post_WL.csv")
head(ax)

ax <- ax[-407,]   #batch 10, typo of ids.

# start times
zt <- read_csv("manuscript/behavior/Post_starttimes.csv")
zts <- zt[zt$times=="startTime",]
ax<- full_join(ax, zts %>% select(batch, starttime = ztime))

ax %>% 
  group_by(batch) %>% 
  filter(date == min(date)) %>% 
  filter(ztime < (4200+starttime)) -> ax1


## Lets start by getting total aggression in each cage
ax1<- ax1 %>% dplyr::select(post_batchcage, batch,Cage,winner,loser,score,group,period,time,ztime)
head(ax1)

##Just getting reorganized cages at 1 hr 
ax1 <- ax1 %>% filter(group == "reorganized") %>% filter(time == "1 hour")
unique(ax1$post_batchcage)
head(ax1)


## split dataframe by individual cage 
l <- split(ax1, ax1$post_batchcage)
lapply(l, head)


## gettings total aggression in each Post Cage
l <- lapply(l, function(x) x) %>% 
  map(~group_by(.,Cage)) %>% 
  map(~arrange(.,ztime)) %>% 
  map(~mutate(.,value = row_number())) %>%
  map(~mutate(.,totAgg = max(value)) )%>% 
  map(~select(.,-value)) %>%
  map(~ungroup(.,))
 

## gettings total aggression in each individual winner 
l2 <- lapply(l, function(x) x) %>% 
  map(~mutate(.,post_id = substrRight(winner,3))) %>% 
  map(~mutate(.,post_id =paste(post_id,post_batchcage))) %>% 
  map(~group_by(.,post_id,winner)) %>%
  map(~arrange(.,ztime)) %>% 
  map(~mutate(.,value = row_number())) %>%
  map(~mutate(.,Agiven = max(value))) %>% 
  map(~ungroup(.,)) %>% 
  map(~select(.,post_batchcage,post_id,time,batch,group,totAgg,Agiven,ztime, -value))
  
lapply(l2, head)

## gettings total recieved aggresssion 
l<- lapply(l, function(x) x) %>% 
  map(~mutate(.,post_id = substrRight(loser,3))) %>% 
  map(~mutate(.,post_id =paste(post_id,post_batchcage))) %>% 
  map(~group_by(.,post_id,loser)) %>% 
  map(~arrange(.,ztime)) %>% 
  map(~mutate(.,value = row_number())) %>%
  map(~mutate(.,ARec = max(value)))%>% 
  map(~ungroup(.,)) %>% 
  map(~select(.,post_batchcage,post_id,time,batch,group,totAgg,ARec, -value)) 


LL <- list(l, l2)
res <- do.call(c, LL)

## add in batch/cage id, and if alpha/beta/gamma/delta cage...
results <- Map(cbind, res, type = c("alphas", "betas", "gammas", "deltas")[as.numeric(substrRight(names(l),1))])
lapply(results, head)


## Make one dataframe
DF <-results %>% 
  map_df(as_tibble)
head(DF)


agg <- DF %>% select(post_batchcage,post_id,time,batch,group,totAgg,Agiven,ARec,type) %>% 
  unique(.)

head(agg)
unique(agg$post_id) # 122 correct 


agg$post_idbatch <- substr(agg$post_id, 1,12)
agg$post_idbatch <- gsub(" ", "", agg$post_idbatch)
agg1 <- agg %>% filter(type != "betas") %>%  select(post_batchcage,post_idbatch,post_id,batch,group,totAgg,Agiven,ARec,type)
head(agg1)

agg1 <- agg1 %>%   pivot_longer(cols = 7:8, names_to="Agg")  

x <-agg1 %>% filter(value!= "NA")
unique(x$post_idbatch) #84 good got rid of betas


dc <- read_csv("manuscript/cort/FullCort.csv")
head(dc)
colnames(dc)


dcx <- dc %>% filter(group == "reorganized") %>% 
  filter(period == "Post") %>% 
  filter(time == "1 hr") %>% 
  select(pre_idbatchcage, post_idbatch,batch, group,period,time,plate,mean_con_ng_ul,Prerank,Postrank,condition) %>% 
  unique(.)


dumb <- x%>% 
  full_join(dcx)


x <-dumb %>% filter(mean_con_ng_ul != "NA")
head(x)
x <- unique(x)

#figuring out werid 1s in Agiven 
x %>% filter(post_id == "3-3 Batch 8 Cage 3")
x %>% filter(post_id == "5-3 Batch 9 Cage 1")

table(x$Postrank)

x$value <- ifelse(x$type == "alphas" & x$Agg == "Agiven" & x$value == 1, "", x$value)
x <-x %>% filter(value != "")


ggplot(x, aes(totAgg,mean_con_ng_ul, color = type))+
         geom_point()

table(x$Postrank)
unique(x$post_idbatch) # 49 great that is what 

head(x)

x <- x %>%   pivot_wider(values_from = value, names_from = Agg)


str(x)
x[is.na(x)] <-"0"
x$Agiven

x <- x %>% group_by(post_idbatch,Agiven,ARec) %>% 
  mutate(.,tot.ind = as.numeric(Agiven) + as.numeric(ARec))

write.csv(x, "manuscript/cort/cortAgg70min.csv", row.names = F)

ggplot(x, aes(as.numeric(totAgg), mean_con_ng_ul, color = as.factor(type)))+
  geom_point(size = 2)+
  labs(title = "",
       x = "Total Aggression",
       y = "Corticosterone (ng/ml)",
       color = "Post Cage") +
  scale_color_manual(values = viridis::viridis(3)) +
  newggtheme


ggplot(x, aes(as.numeric(Agiven), mean_con_ng_ul, color = as.factor(Postrank)))+
  geom_point()


ggplot(x, aes(as.numeric(ARec), mean_con_ng_ul, color = as.factor(Postrank)))+
  geom_point()
  

x$type <- factor(x$type, levels = c("alphas", "gammas", "deltas"))
x$domgroup <- ifelse(x$Postrank == 1, "Dominant", "Subordinate")

ggplot(x, aes(tot.ind, mean_con_ng_ul, color = as.factor(domgroup)))+
  geom_point(size = 3, alpha = .6)+
  # geom_smooth(method= "lm",se=F)+
  labs(title = "",
       x = "Total Individual Aggression",
       y = "Corticosterone (ng/ml)",
       color = "Post Domgroup") +
  scale_color_manual(values = viridis::viridis(3)) +
  facet_wrap(~type)+
theme_test()



## linear_models
library(MASS)
library(lme4)
library(emmeans)
library(car)
colnames(x)

x$type <- as.factor(x$type)
x$tot.ind <- as.numeric(x$tot.ind)
x$domgroup <- as.factor(x$domgroup)

agg.m <- glmer(mean_con_ng_ul ~ tot.ind + type  +(1|post_batchcage)+ (1|post_idbatch)+ (1|batch)+(1|plate), data =x, family= Gamma(link ="inverse"))
summary(agg.m)

AIC(agg.m)
AIC(agg.m1)
AIC(agg.m2)



## why are these giving me bery different results -- maybe I need to talk another stat class
agg.m1 <- glmer(mean_con_ng_ul ~ tot.ind + type +domgroup +(1|post_batchcage)+ (1|post_idbatch)+ (1|batch)+(1|plate), data =x, family= Gamma(link ="inverse"))
summary(agg.m1)

agg.m2 <- glmer(mean_con_ng_ul ~ tot.ind +domgroup+type+(1|post_batchcage/batch)+ (1|post_idbatch/batch)+ (1|plate), data =x, family=Gamma(link="log"))
summary(agg.m2)

blah <- emmeans(agg.m2, ~ tot.ind + domgroup + type)
pairs(blah)
blah

##assumptions: 
## all terrible 
acf(resid(agg.m2))  
qqPlot(resid(agg.m2))
hist(resid(agg.m2))

plot(x=fitted(agg.m2), y = resid(agg.m2))

durbinWatsonTest(resid(agg.m2))
shapiro.test(resid(agg.m2))



##### analysis at 25 hr

ax %>% 
  group_by(batch) %>% 
  filter(date == min(date)) %>% 
  filter(ztime < (4200+starttime)) -> ax1


## Lets start by getting total aggression in each cage
ax1<- ax1 %>% dplyr::select(post_batchcage, batch,Cage,winner,loser,score,group,period,time,ztime)
head(ax1)

##Just getting reorganized cages at 1 hr 
ax25 <- ax1 %>% filter(group == "reorganized") %>% filter(time == "25 hours")
unique(ax25$post_batchcage)
head(ax25)


## split dataframe by individual cage 
l <- split(ax25, ax25$post_batchcage)
lapply(l, head)


## gettings total aggression in each Post Cage
l <- lapply(l, function(x) x) %>% 
  map(~group_by(.,Cage)) %>% 
  map(~arrange(.,ztime)) %>% 
  map(~mutate(.,value = row_number())) %>%
  map(~mutate(.,totAgg = max(value)) )%>% 
  map(~select(.,-value)) %>%
  map(~ungroup(.,))


## gettings total aggression in each individual winner 
l2 <- lapply(l, function(x) x) %>% 
  map(~mutate(.,post_id = substrRight(winner,3))) %>% 
  map(~mutate(.,post_id =paste(post_id,post_batchcage))) %>% 
  map(~group_by(.,post_id,winner)) %>%
  map(~arrange(.,ztime)) %>% 
  map(~mutate(.,value = row_number())) %>%
  map(~mutate(.,Agiven = max(value))) %>% 
  map(~ungroup(.,)) %>% 
  map(~select(.,post_batchcage,post_id,time,batch,group,totAgg,Agiven,ztime, -value))

lapply(l2, head)

## gettings total recieved aggresssion 
l<- lapply(l, function(x) x) %>% 
  map(~mutate(.,post_id = substrRight(loser,3))) %>% 
  map(~mutate(.,post_id =paste(post_id,post_batchcage))) %>% 
  map(~group_by(.,post_id,loser)) %>% 
  map(~arrange(.,ztime)) %>% 
  map(~mutate(.,value = row_number())) %>%
  map(~mutate(.,ARec = max(value)))%>% 
  map(~ungroup(.,)) %>% 
  map(~select(.,post_batchcage,post_id,time,batch,group,totAgg,ARec, -value)) 


LL <- list(l, l2)
res <- do.call(c, LL)

## add in batch/cage id, and if alpha/beta/gamma/delta cage...
results <- Map(cbind, res, type = c("alphas", "betas", "gammas", "deltas")[as.numeric(substrRight(names(l),1))])
lapply(results, head)


## Make one dataframe
DF <-results %>% 
  map_df(as_tibble)
head(DF)


agg <- DF %>% select(post_batchcage,post_id,time,batch,group,totAgg,Agiven,ARec,type) %>% 
  unique(.)

head(agg)
unique(agg$post_id) # 96 correct 


agg$post_idbatch <- substr(agg$post_id, 1,12)
agg$post_idbatch <- gsub(" ", "", agg$post_idbatch)
agg1 <- agg %>% filter(type != "betas") %>%  select(post_batchcage,post_idbatch,post_id,batch,group,totAgg,Agiven,ARec,type)
head(agg1)

agg1 <- agg1 %>%   pivot_longer(cols = 7:8, names_to="Agg")  

x <-agg1 %>% filter(value!= "NA")
unique(x$post_idbatch) #70, but should be 72 -_- missing 2 animals after getting rid of betas
head(x)

dc <- read_csv("manuscript/cort/FullCort.csv")
head(dc)
colnames(dc)


cort <- dc %>% filter(group == "reorganized") %>% 
  filter(period == "Post") %>% 
  filter(time == "25hr") %>% 
  select(pre_idbatchcage, post_idbatch,batch, group,period,time,plate,mean_con_ng_ul,Prerank,Postrank,condition) %>% 
  unique(.)
head(cort)

day1 <- x %>% 
  full_join(cort)
head(day1)

x25 <-day1 %>% filter(mean_con_ng_ul != "NA")
head(x25)

head(x25)
table(x25$Postrank)

x25$value <- ifelse(x25$type == "alphas" & x25$Agg == "Agiven" & x25$value == 1, "", x25$value)
x25 <-x25 %>% filter(value != "")


ggplot(x25, aes(totAgg,mean_con_ng_ul, color = type))+
  geom_point()

table(x$Postrank)
unique(x$post_idbatch) # 70, still missing 2 

head(x)

x25 <- x25 %>%   pivot_wider(values_from = value, names_from = Agg)


str(x25)
x25[is.na(x25)] <-"0"
x25$Agiven

x25 <- x25%>% group_by(post_idbatch,Agiven,ARec) %>% 
  mutate(.,tot.ind = as.numeric(Agiven) + as.numeric(ARec))

write.csv(x, "manuscript/cort/cortAgg25hr.csv", row.names = F)


ggplot(x25, aes(as.numeric(totAgg), mean_con_ng_ul, color = as.factor(type)))+
  geom_point(size = 2)+
  labs(title = "",
       x = "Total Aggression",
       y = "Corticosterone (ng/ml)",
       color = "Post Cage") +
  scale_color_manual(values = viridis::viridis(3))+
 theme_classic()


ggplot(x25, aes(as.numeric(Agiven), mean_con_ng_ul, color = as.factor(Postrank)))+
  geom_point()+
  theme_classic()


ggplot(x25, aes(as.numeric(ARec), mean_con_ng_ul, color = as.factor(Postrank)))+
  geom_point()+
  theme_classic()


x25$type <- factor(x25$type, levels = c("alphas", "gammas", "deltas"))
x25$domgroup <- ifelse(x25$Postrank == 1, "Dominant", "Subordinate")

ggplot(x25, aes(tot.ind, mean_con_ng_ul, color = as.factor(domgroup)))+
  geom_point(size = 3, alpha = .6)+
  # geom_smooth(method= "lm",se=F)+
  labs(title = "",
       x = "Total Individual Aggression",
       y = "Corticosterone (ng/ml)",
       color = "Post Domgroup") +
  scale_color_manual(values = viridis::viridis(3)) +
  facet_wrap(~type)+
  theme_test()

