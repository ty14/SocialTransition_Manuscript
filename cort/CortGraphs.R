# Graphs for the paper as organized as possible. 


#libraries and source
library(tidyverse)
library(viridis)

df <- read_csv("manuscript/cort/FullCort.csv")
head(df)

df <- df %>% select(-Well, -final_conc_ng_ul)

## just getting one data point 
df1 <- unique(df) 

# first graph boxplot of CORT (yaxis)  vs  group (to-be-reorganized  vs control on x-axis)
g1 <- df %>%filter(period == "Pre")


g1$group <- ifelse(g1$group == "reorganized", "To-Be-Reorganized","Control")


ggplot(g1, aes(group, mean_con_ng_ul,color =group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21,
                 alpha = 0.8,
                 jitter.height = 0.02, jitter.width = 0.030, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  labs(title = "",
       x = "",
       y = "Corticosterone (ng/ml)") +
  scale_color_manual(values = viridis::viridis(4)) +
  scale_fill_manual(values = viridis::viridis(4))+
  newggtheme+
  theme( legend.position = 'none')


## descriptive stats 
g1 %>% group_by(group) %>%
  summarise(median = median(mean_con_ng_ul),
            lqr = quantile(mean_con_ng_ul,.25),
            uqr = quantile(mean_con_ng_ul,.75))


## second graph no rank differences 
g1$Prerank <- ifelse(g1$Prerank == 1, "Alphas", g1$Prerank)
g1$Prerank <- ifelse(g1$Prerank == 3,"Gammas", g1$Prerank)
g1$Prerank <- ifelse(g1$Prerank == 4, "Deltas", g1$Prerank)
g1$Prerank <-factor(g1$Prerank, levels = c("Alphas","Gammas", "Deltas"))

ggplot(g1, aes(Prerank, mean_con_ng_ul,color =Prerank))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21,
                 alpha = 0.8,
                 jitter.height = 0.02, jitter.width = 0.030, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  labs(title = "",
       x = "",
       y = "Corticosterone (ng/ml)") +
  scale_color_manual(values = viridis::viridis(4)) +
  scale_fill_manual(values = viridis::viridis(4))+
  newggtheme+
  theme( legend.position = 'none')


## descriptive stats for prerank  
g1 %>% group_by(Prerank) %>%
  summarise(median = median(mean_con_ng_ul),
            lqr = quantile(mean_con_ng_ul,.25),
            uqr = quantile(mean_con_ng_ul,.75))



### Pre and Post data
head(df1)
df1$time <- ifelse(df1$time == "1 hr", "Reorganized 70 min", df1$time)
df1$time <- ifelse(df1$time == "25hr", "Reorganized 25 hr", df1$time)
## Time for the control data 
df1$time <- ifelse(df1$time == "Reorganized 70 min" & df1$group == "control","Control 70 min", df1$time)
df1$time <- ifelse(df1$time == "Reorganized 25 hr" & df1$group == "control","Control 25 hr", df1$time)

#Making things factors
df1$period <- factor(df1$period, levels = c("Pre", "Post"))
df1$time <- factor(df1$time,levels =c ("Reorganized 70 min", "Control 70 min", "Reorganized 25 hr", "Control 25 hr"))

## Creating domgroup based on Prerank 
df1$domgroup <- ifelse(df1$Prerank == 3|df1$Prerank==4, "Subordinate", "Dominant") 


ggplot(df1, aes(period, mean_con_ng_ul))+
  geom_point(aes(color =time), alpha = .6, size = 3) +
  geom_line(aes(group = pre_idbatchcage, color =time), size = .5) +
  labs(x = "",
       y = "Corticosterone (ng/ml)") +
  scale_color_manual(values = viridis::viridis(4)) +
  facet_grid(time~domgroup, space="free", scales="free")+
  newggtheme 


## getting things
post <- df1 %>% 
  group_by(domgroup,time,period,group) %>%
  summarize(mean_post = mean(mean_con_ng_ul),
            sd_post = sd(mean_con_ng_ul),
            n= n(),
            median_post = median(mean_con_ng_ul)) %>%
  mutate(semx = sd_post/sqrt(n)) %>% 
  filter(!is.na(semx)) %>% 
  mutate( lower_meanp = mean_post + qt((1-0.95)/2, n - 1) * semx,
          upper_meanp = mean_post - qt((1-0.95)/2, n - 1) * semx) %>%
  mutate( lower_medp = median_post + qt((1-0.95)/2, n - 1) * semx,
          upper_medp = median_post - qt((1-0.95)/2, n - 1) * semx) %>%
  ungroup()


post <-df1 %>% dplyr::select(pre_idbatchcage,domgroup,period, time,mean_con_ng_ul) %>% full_join(post)


## Real graph 
ggplot(post, aes(x=period , y = mean_con_ng_ul)) +
  geom_ribbon(aes(ymin = lower_medp,
                  ymax = upper_medp, group=domgroup,fill=domgroup,alpha=5)) +
  geom_line(aes(y=median_post, group = domgroup, color = domgroup), size=1.5)+
  geom_line(aes(group = pre_idbatchcage), color= "gray", size = .6, alpha=.5) +
  geom_point(aes(group =domgroup), color ="gray", alpha = .5, size = .6)+
  ylim(0, 600)+
  labs(   x = "",
          y=  "Corticosterone (ng/ml)") +
  scale_color_manual(values = viridis::viridis(3)) +
  facet_grid(time~domgroup, space="free", scales="free")+
  newggtheme +
  theme(legend.position = "none")

#descriptive stats
post %>% group_by(domgroup,time) %>%
  summarise(median = median(mean_con_ng_ul),
            lqr = quantile(mean_con_ng_ul,.25),
            uqr = quantile(mean_con_ng_ul,.75))

## Looking at the difference between pre and post data 


## looking at how condition maps 
head(df)

df1 <- df1 %>% filter(time == "Reorganized 70 min") 

df1$condition1<- ifelse(df1$condition == "ascenders",  "SUB -> DOM", df1$condition)
df1$condition1 <-ifelse(df1$condition == "descenders",  "DOM -> SUB", df1$condition1)
df1$condition1 <- ifelse(df1$condition == "same" & df1$domgroup== "Dominant", "DOM -> DOM", df1$condition1)
df1$condition1 <- ifelse(df1$condition == "same" & df1$domgroup== "Subordinate", "SUB -> SUB", df1$condition1)
df1$condition1 <- factor(df1$condition1, levels = c("DOM -> DOM","DOM -> SUB","SUB -> SUB","SUB -> DOM"))
sum <-  df1 %>% 
  group_by(domgroup,time,period,condition1) %>%
  summarize(sdx = sd(mean_con_ng_ul),
            n= n(),
            medianx = median(mean_con_ng_ul)) %>%
  mutate(semx = sdx/sqrt(n)) %>% 
  filter(!is.na(semx)) %>%
  mutate( lower_med = medianx + qt((1-0.95)/2, n - 1) * semx,
          upper_med = medianx - qt((1-0.95)/2, n - 1) * semx) %>%
  ungroup()

one <- df1 %>% dplyr::select(pre_idbatchcage,domgroup,period, time,mean_con_ng_ul,condition, condition1)

one <- one %>%full_join(sum)

xx <-one %>% select(pre_idbatchcage,mean_con_ng_ul, period,condition1, domgroup) %>% unique(.) %>% 
  group_by(pre_idbatchcage) %>% 
  pivot_wider(values_from = mean_con_ng_ul, names_from = period) %>% 
  mutate(diff = Post - Pre) %>% ungroup()


xx %>% group_by(condition1) %>%
  summarise(median = median(diff),
            lqr = quantile(diff,.25),
            uqr = quantile(diff,.75))


ggplot(one, aes(period,mean_con_ng_ul))+
  geom_ribbon(aes(ymin = lower_med,
                  ymax = upper_med, group=domgroup,fill=domgroup, alpha = 5)) +
  scale_alpha(guide = 'none')+
  scale_fill_discrete(guide=FALSE)+
  geom_line(aes(group = pre_idbatchcage),color ="gray", size = .4, alpha=.4) +
  geom_point(aes(group =domgroup),color="gray", alpha = .4, size = 1)+
  geom_line(aes(y=medianx, group = domgroup, color = domgroup), size=1)+
  facet_wrap(~condition1)+
  labs(x = "",
       y=  "Corticosterone (ng/ml)",
       color = "Pre Domgroup") +
  scale_color_manual(values = viridis::viridis(3))+
  newggtheme+
  theme(legend.position = "none")



# Rank differences at 1 hr post -  between cages (alpha/beta/gamma/delta)  
# and between ranks (i.e. 1s/3s/4s)

post1<- read_csv("manuscript/cort/cortAgg.csv")
head(post1)

post1 <- post %>% dplyr::select(-totAgg,-tot.ind,-Agiven, -ARec, -post_id)
head(post1)

ggplot(df)



# Rank differences at 25 hr post -  between cages (alpha/beta/gamma/delta)  
# and between ranks (i.e. 1s/3s/4s)
head(df)
df25$Prerank

df25 <-df %>% filter(group == "reorganized")
df25$type<- ifelse(df25$Prerank == 1, "alphas", "")
df25$type<- ifelse(df25$Prerank == 3, "gammas", df25$type)
df25$type<- ifelse(df25$Prerank == 4, "deltas", df25$type)

df25$type <-factor(df25$type, levels = c("alphas", "gammas", "deltas"))
df25$time <- ifelse(df25$time == "1 hr", "70 min", df25$time)
df25$time <- ifelse(df25$time == "25hr", "25 hr", df25$time)
df25$time <- factor(df25$time, levels = c("70 min", "25 hr"))

df25 %>% filter(period == "Post") %>% 
ggplot(., aes(type,mean_con_ng_ul, color =as.factor(Postrank)))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21,
                 alpha = 0.8,
                 jitter.height = 0.02, jitter.width = 0.030, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  labs(title = "",
       x = "New Cage",
       y = "Corticosterone (ng/ml)",
       color = "Post Rank") +
  scale_color_manual(values = viridis::viridis(4)) +
  scale_fill_manual(values = viridis::viridis(4))+
  facet_wrap(~time)+
  newggtheme



library(lme4)
library(emmeans)
colnames(df25)

df25$Postrank <- as.factor(df25$Postrank)
df25$domgroupPost <- as.factor(df25$domgroupPost)

df25x <- df25 %>% filter(time == "25 hr")
df1 <- df25 %>%  filter(time == "70 min")


colnames(df25)
glm <- glmer(mean_con_ng_ul ~ domgroupPost + type+ (1|pre_idbatchcage) + (1|pre_id)+ (1|post_id) +(1|batch) + (1|plate), data = df25x, family = Gamma(link= "log"))
summary(glm)

x <- emmeans(glm, ~ domgroupPost + type)
pairs(x)


glm1 <- glmer(mean_con_ng_ul ~ type + domgroupPost+ (1|pre_idbatchcage) + (1|pre_id)+ (1|post_id) +(1|batch) + (1|plate), data = df1, family = Gamma(link= "log"))
summary(glm1)

x1 <- emmeans(glm1, ~ domgroupPost + type)
pairs(x1)


df1 %>% group_by(type) %>% 
  summarise(median = median(mean_con_ng_ul),
            lqr = quantile(mean_con_ng_ul,.25),
            uqr = quantile(mean_con_ng_ul,.75))


df25x %>% group_by(type, domgroupPost) %>% 
  summarise(median = median(mean_con_ng_ul),
            lqr = quantile(mean_con_ng_ul,.25),
            uqr = quantile(mean_con_ng_ul,.75))


### random thing 
df$domgroup <- ifelse(df$Prerank == 1, "Dominant","Subordinant")
df$domgroupPost <- ifelse(df$Postrank == 1, "Dominant","Subordinant")
head(df)

df <- df %>% dplyr::select(-Well,-dilution_factor,-final_conc_ng_ul) %>%  unique(.)

domsub <- df %>%filter(domgroup == "Dominant" & domgroupPost == "Subordinant") #51%
dom <- df %>%filter(domgroup == "Dominant" & domgroupPost == "Dominant")#49%
dom1 <- df %>%filter(domgroup == "Dominant")


subdom <- df %>%filter(domgroupPost == "Dominant" & domgroup == "Subordinant")#51 35 %
sub <- df %>%filter(domgroup == "Subordinant"& domgroupPost == "Subordinant")#98 #65 %
sub1 <- df %>%filter(domgroup == "Subordinant")#149

