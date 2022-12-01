# cort analysis 
library(viridis)
library(glue)
library(lme4)
library(MASS)
library(car)
library(tidyverse)

#source
df <- read_csv("manuscript/cort/FullCort.csv") 
head(df)

df<- df %>% dplyr::select(-Well, -final_conc_ng_ul) %>%unique(.)


#changing times
df$time <- ifelse(df$time == "1 hr", "70 min", df$time)
df$time <- ifelse(df$time == "25hr", "25 hr", df$time)
table(df$time)

##Changing period
head(df)
df$period <- ifelse(df$period == "Pre", 0, df$period)
df$period <- ifelse(df$period == "Post" & df$time == "70 min", 1, df$period)
df$period <- ifelse(df$period == "Post" & df$time == "25 hr", 2, df$period)

df <- df %>% 
  dplyr::select( -time)

##adding domgroup 
df$domgroup <- ifelse(df$Prerank == 3|df$Prerank==4, "Subordinate", "Dominant")  
## post domgroup 
df$domgroupPost <- ifelse(df$Postrank == 3|df$Postrank==4, "Subordinate", "Dominant")  

## Getting rid of duplicates 
df <- unique(df)


## making sure everything is a factor ( idk if you need to do this, but I think so)
df$period <- as.factor(df$period)
df$Prerank <- factor(df$Prerank, levels = c("1","3", "4"))
df$Postrank <- factor(df$Postrank, levels = c("1", "3","4"))
df$batch <- as.factor(df$batch)
df$plate <- as.factor(df$plate)
df$pre_idbatchcage <- as.factor(df$pre_idbatchcage)
df$condition<- factor(df$condition, levels = c("control", "same", "ascenders", "descenders"))
df$post_idbatch <- as.factor(df$post_idbatch)
df$group <- as.factor(df$group)


# Figure 1
g1 <- df %>%filter(period == "0")
g1$group <- ifelse(g1$group == "reorganized", "Reorganized","Control")
head(g1)

p1 <- ggplot(g1, aes(group, mean_con_ng_ul,color =group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21,
                 alpha = 0.4,
                 jitter.height = 0.02, jitter.width = 0.030, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  labs(title = "",
       x = "",
       y = "Corticosterone (ng/ml)") +
  scale_color_manual(values = viridis::viridis(4)) +
  scale_fill_manual(values = viridis::viridis(4)) +theme_classic() +
  theme(axis.text.x = element_text(vjust = 1, size = 15),
        legend.position = 'none',
        axis.text.y = element_text(hjust = 0.5, size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title = element_text(size = 20),
        axis.text= element_text(color="#3C3C3C", size=20),
        strip.background = element_blank() 
  )



#Figure 2 

g1$Prerank <- ifelse(g1$Prerank == 1, "Alphas", g1$Prerank)
g1$Prerank <- ifelse(g1$Prerank == 2,"Gammas", g1$Prerank)
g1$Prerank <- ifelse(g1$Prerank == 3, "Deltas", g1$Prerank)
g1$Prerank <-factor(g1$Prerank, levels = c("Alphas","Gammas", "Deltas"))

p2 <-ggplot(g1, aes(Prerank, mean_con_ng_ul,color =Prerank, fill = Prerank))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21,
                 alpha = 0.4,
                 jitter.height = 0.02, jitter.width = 0.030, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  labs(title = "",
       x = "",
       y = "Corticosterone (ng/ml)") +
  scale_color_manual(values = viridis::viridis(4)) +
  scale_fill_manual(values = viridis::viridis(4))+theme_classic() +
  theme(axis.text.x = element_text(vjust = 1, size = 15),
        legend.position = 'none',
        axis.text.y = element_text(hjust = 0.5, size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title = element_text(size = 20),
        axis.text= element_text(color="#3C3C3C", size=20),
        strip.background = element_blank() 
  )


pm3 <- grid.arrange(p1,p2,nrow=1)
ggsave("manuscript/behavior/results_figures/Cort_pre.png",pm3,height =5, width =10, dpi=600)

##Pre reorganization stats for Figure 1 and 2 

hist(g1$mean_con_ng_ul) #not normal 
g1$Prerank

ggplot(g1, aes(Prerank,mean_con_ng_ul, color = Prerank)) +
 geom_boxplot()+
  geom_jitter() +
  facet_wrap(~batch)

## model for the 2 pre graphs - no differences obivously 
colnames(g1)

g1$Prerank <- factor(g1$Prerank, levels = c("Gammas","Deltas", "Alphas"))
pre.glm <-glmer(mean_con_ng_ul~group+ Prerank+(1|batch)+(1|plate)+(1|pre_id), data =g1,family = Gamma(link = "log"))
summary(pre.glm)
AIC(pre.glm)


#checks 
acf(resid(pre.glm))
qqPlot(resid(pre.glm))
hist(resid(pre.glm))
plot(pre.glm)

durbinWatsonTest(resid(pre.glm))
shapiro.test(resid(pre.glm))


#Figure 3 comparing pre and post cort  by pre-domgroup 
df1 <- read_csv("manuscript/cort/FullCort.csv") 
head(df1)
df1<- df1 %>% dplyr::select(-Well, -final_conc_ng_ul) %>%unique(.)
df1$time <- ifelse(df1$time == "1 hr", "Reorganized 70 min", df1$time)
df1$time <- ifelse(df1$time == "25hr", "Reorganized 25 hr", df1$time)
## Time for the control data 
df1$time <- ifelse(df1$time == "Reorganized 70 min" & df1$group == "control","Control", df1$time)
df1$time <- ifelse(df1$time == "Reorganized 25 hr" & df1$group == "control","Control", df1$time)

#Making things factors
df1$period <- factor(df1$period, levels = c("Pre", "Post"))
df1$time <- factor(df1$time,levels =c ("Reorganized 70 min", "Reorganized 25 hr", "Control"))
##adding domgroup 
df1$domgroup <- ifelse(df1$Prerank == 3|df1$Prerank==4, "Subordinate", "Dominant")  

## getting stats
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


## Real Figure 
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
  # newggtheme +
  theme(legend.position = "none", 
        strip.text = element_text(size = 10))


  cort4 <- ggplot(post, aes(x=period , y = mean_con_ng_ul)) +
  # geom_ribbon(aes(ymin = lower_medp,
  #                 ymax = upper_medp, group=domgroup,fill=domgroup,alpha=5)) +
 # geom_line(aes(y=median_post, group = domgroup, color = domgroup), size=1.5)+
  geom_line(aes(group = pre_idbatchcage, color=domgroup), size = .75, alpha = .75) +
  geom_point(aes(color=domgroup), alpha = .75, size = 1)+
  scale_color_manual(values = c("#238A8DFF", "#FDE725FF"), name = "PreRank",  labels=c("DOM","SUB"))+
  ylim(0, 600)+
  labs(   x = "",
          y=  "Corticosterone (ng/ml)") +
  facet_wrap(~time,nrow =1)+
  theme(axis.text.x = element_text(vjust = 1, size = 15),
        axis.text.y = element_text(hjust = 0.5, size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title = element_text(size = 20),
        axis.text= element_text(color="#3C3C3C", size=15),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15)
  )


ggsave("manuscript/behavior/results_figures/Cort_7025con.png",cort4,height =4, width =8, dpi=600)

# scale_fill_manual(values = c("#404788FF", "#238A8DFF", "#55C667FF", "#FDE725FF"),
#                  name="Rank", labels=c("1","2","3","4")) +
#   scale_color_manual(values = c("#404788FF", "#238A8DFF", "#55C667FF", "#FDE725FF"),
#                      name="Rank", labels=c("1","2","3","4"))

## mix model for Figure 3 
colnames(df1)
df1$domgroup <- as.factor(df1$domgroup)
df1$time2 <- factor(df1$time, levels = c( "Control 70 min","Reorganized 70 min","Control 25 hr", "Reorganized 25 hr"))
prepost<- glmer(mean_con_ng_ul~period + time+domgroup +(1|pre_idbatchcage)+ (1|post_idbatch)+(1|batch)+(1|plate), data =df1, family = Gamma(link = "log"))
summary(prepost )
AIC(prepost) #better 
broom::tidy(prepost)
prepost1<- glmer(mean_con_ng_ul~period + time2+domgroup +(1|pre_idbatchcage)+ (1|post_idbatch)+(1|batch)+(1|plate), data =df1, family = Gamma(link = "log"))
summary(prepost1)

colnames(df1)

r <- df1 %>% filter(time2 == "Reorganized 70 min")
pp_reorg<- glmer(mean_con_ng_ul~period +(1|pre_idbatchcage)+ (1|post_idbatch)+(1|batch)+(1|plate), data = r, family = Gamma(link = "log"))
summary(pp_reorg)

c <- df1 %>% filter(time2 == "Control 70 min")
pp_c<- glmer(mean_con_ng_ul~period +(1|pre_idbatchcage)+ (1|post_idbatch)+(1|batch)+(1|plate), data = c, family = Gamma(link = "log"))
summary(pp_c)

r25 <- df1 %>% filter(time2 != "Reorganized 70 min") %>% filter(time2 != "Control 70 min")
pp_25<- glmer(mean_con_ng_ul~time2 +(1|pre_idbatchcage)+ (1|post_idbatch)+(1|batch)+(1|plate), data = r25, family = Gamma(link = "log"))
summary(pp_25)

#check
acf(resid(prepost))
qqPlot(resid(prepost))
hist(resid(prepost))
plot(prepost)

durbinWatsonTest(resid(prepost))
shapiro.test(resid(prepost))


# or look at the diff in mixed model 
dp <-df1 %>% 
  group_by(pre_idbatchcage) %>% 
  pivot_wider(values_from = mean_con_ng_ul, names_from = period) %>% 
  mutate(diff = Post - Pre) %>% ungroup() %>% 
  pivot_longer(cols = 19:20, names_to="period")
head(dp)
# dp$diff <- as.numeric(dp$diff)


#mixed model 2 this won't let me add in domgroup 
# hist(dp$diff)
dp$time
dp1 <- df1 %>% filter(time == "Reorganized 25 hr") 
# dp1 <- dp1 %>% dplyr::select(domgroup, condition,time,pre_idbatchcage,batch,plate)
dp1$condition

dp <- df1 %>% filter(time == "Reorganized 70 min") 

dp$condition1<- ifelse(dp$condition == "ascenders",  "Ascenders", dp$condition)
dp$condition1 <-ifelse(dp$condition == "descenders",  "Descenders", dp$condition1)
dp$condition1 <- ifelse(dp$condition == "same" & dp$domgroup== "Dominant", "Dominant", dp$condition1)
dp$condition1 <- ifelse(dp$condition == "same" & dp$domgroup== "Subordinate", "Subordinate", dp$condition1)
dp$condition1 <- factor(dp$condition1, levels = c("Descenders","Dominant","Subordinate","Ascenders"))

dp1$condition1<- ifelse(dp1$condition == "ascenders",  "Ascenders", dp1$condition)
dp1$condition1 <-ifelse(dp1$condition == "descenders",  "Descenders", dp1$condition1)
dp1$condition1 <- ifelse(dp1$condition == "same" & dp1$domgroup== "Dominant", "Dominant", dp1$condition1)
dp1$condition1 <- ifelse(dp1$condition == "same" & dp1$domgroup== "Subordinate", "Subordinate", dp1$condition1)
dp1$condition1 <- factor(dp1$condition1, levels = c("Descenders","Dominant","Subordinate","Ascenders"))



colnames(dp)

ggplot(dp, aes(condition1, diff)) +
  geom_boxplot()+
  geom_jitter()

diff70 <- ggplot(dp, aes(condition1, mean_con_ng_ul,color =condition1, fill = condition1))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21,
                 alpha = 0.4,
                 jitter.height = 0.02, jitter.width = 0.030, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  labs(title = "70 min",
       x = "",
       y = "Post Corticosterone (ng/ml)") +
  scale_y_continuous(limits = c(0,550), breaks= c(0,100,300,500))+
  scale_color_manual(values = viridis::viridis(4)) +
  scale_fill_manual(values = viridis::viridis(4))+
  theme(axis.text.x = element_text(vjust = 1, size = 15),
        axis.text.y = element_text(hjust = 0.5, size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = 'none',
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title = element_text(size = 20),
        axis.text= element_text(color="#3C3C3C", size=15),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        element_text(size = 20)
  )
  



diff25 <- ggplot(dp1, aes(condition1, mean_con_ng_ul,color =condition1, fill = condition1))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21,
                 alpha = 0.4,
                 jitter.height = 0.02, jitter.width = 0.030, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  labs(title = "25 hr",
       x = "",
       y = "") +
  scale_y_continuous(limits = c(0,550), breaks= c(0,100,300,500))+
  scale_color_manual(values = viridis::viridis(4)) +
  scale_fill_manual(values = viridis::viridis(4))+
  theme(axis.text.x = element_text(vjust = 1, size = 15),
        axis.text.y = element_text(hjust = 0.5, size = 15),
        panel.grid.major = element_blank(),
        legend.position = 'none',
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title = element_text(size = 20),
        axis.text= element_text(color="#3C3C3C", size=15),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        element_text(size = 20)
  )

postcort <- gridExtra::grid.arrange(diff70, diff25, nrow =1)
ggsave("manuscript/behavior/results_figures/postcort.png",postcort,height =5, width =11, dpi=600)

hist(dp1$diff)
dp1 <- na.omit(dp1$diff)

diff.lm<- lmer(diff ~ condition1 + (1|pre_idbatchcage)+(1|batch)+(1|plate), data =dp)
summary(diff.lm)
AIC(diff.lm) #better 
broom::tidy(diff.lm)



dp1$condition1<- ifelse(dp1$condition == "Ascenders",1, dp1$condition1)
dp1$condition1 <-ifelse(dp1$condition == "Descenders",2, dp1$condition1)
dp1$condition1 <- ifelse(dp1$condition == "Dominant",3, dp1$condition1)
dp1$condition1 <- ifelse(dp1$condition == "Subordinate", 4,dp1$condition1)

diff.lm25<- lmer(diff ~ condition1 + (1|pre_idbatchcage)+(1|batch)+(1|plate), data =dp1)
summary(diff.lm25)
dp2<- dp1 %>% filter(condition1 != "3") %>% filter(condition1 != "4")
 
t.test (dp2$diff)
#get pvlaues 
# extract coefficients
coefs <- data.frame(coef(summary(diff.lm)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

# Figure 4 -- looking at different conditions at 1hr 
head(df)

dp1 <- df1 %>% filter(time == "Reorganized 70 min") 

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

xx <-one %>% dplyr::select(pre_idbatchcage,mean_con_ng_ul, period,condition1, domgroup) %>% unique(.) %>% 
  group_by(pre_idbatchcage) %>% 
  pivot_wider(values_from = mean_con_ng_ul, names_from = period) %>% 
  mutate(diff = Post - Pre) %>% ungroup()
head(xx)
xx <- xx %>% pivot_longer(cols = 4:5, names_to="period")

xx <- xx %>% select(-period)

xx <- xx %>% full_join(df1)

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

colnames(df1)
diff <- df1 %>% dplyr::select(pre_id,period,mean_con_ng_ul,condition1) 

diff%>% pivot_wider(., names_from = period, values_from = mean_con_ng_ul)
one %>% group_by(condition1)%>% pivot_wider(., names_from = period, values_from = mean_con_ng_ul)
## mix model for Figure 3 
hist(xx$diff)
colnames(xx) 
table(xx$condition1, xx$domgroup)
#this works
diff.lm<- lm(diff~condition1, data =xx)
summary(diff.lm)

cond <- xx
cond$condition1 <- factor(cond$condition1, levels = c("SUB -> DOM", "SUB -> SUB","DOM -> SUB", "DOM -> DOM"))

diff.lm2<- lm(diff~condition1, data =cond)
summary(diff.lm2)

cond2 <- xx
cond2$condition2 <- factor(cond$condition1, levels = c("SUB -> SUB","SUB -> DOM","DOM -> SUB", "DOM -> DOM"))
diff.lm3<- lm(diff~condition2, data =cond2)
summary(diff.lm3)

#this doesn't work -_-
diff.lmer<- lmer(diff~ condition1 +(1|pre_idbatchcage)+(1|post_idbatch)+(1|batch)+(1|plate), data =cond)
summary(diff.lmer)

#this works tho
df1$condition1
df2 <- df1 
df2$condition2 <- factor(df1$condition1, levels = c("SUB -> DOM", "SUB -> SUB","DOM -> SUB", "DOM -> DOM"))

con.glm <- glmer(mean_con_ng_ul~period+ condition1 +(1|pre_idbatchcage)+(1|post_idbatch)+(1|batch)+(1|plate), data =df1, family =Gamma(link = "log"))
summary(con.glm)

con.glm2 <- glmer(mean_con_ng_ul~period+ condition2 +(1|pre_idbatchcage)+(1|post_idbatch)+(1|batch)+(1|plate), data =df2, family =Gamma(link = "log"))
summary(con.glm2)
AIC(con.glm2)

plot(con.glm2)
acf(resid(con.glm2))
qqPlot(resid(con.glm2))

durbinWatsonTest(resid(con.glm2))
shapiro.test(resid(con.glm2))


### other models 

cond <- unique(xx)
colnames(cond)
cond$condition1 <- factor(cond$condition1, levels = c("SUB -> DOM", "SUB -> SUB","DOM -> SUB", "DOM -> DOM"))

diff.lm2<- glmer(mean_con_ng_ul~period+condition1+ (1|pre_idbatchcage)+ (1|post_idbatch)+(1|batch)+(1|plate), data =cond,family = Gamma(link = "log"))
summary(diff.lm2)
plot(diff.lm2)
acf(resid(diff.lm2))
qqPlot(resid(diff.lm2))


cond2 <- unique(xx)
cond2$condition2 <- factor(cond$condition1, levels = c("SUB -> SUB","SUB -> DOM","DOM -> SUB", "DOM -> DOM"))
diff.lm3<- glmer(mean_con_ng_ul~period+condition2+ (1|pre_idbatchcage)+ (1|post_idbatch)+(1|batch)+(1|plate), data =cond2,family = Gamma(link = "log"))
summary(diff.lm3)
plot(diff.lm3)
acf(resid(diff.lm3))
qqPlot(resid(diff.lm3))

cond3 <- unique(xx)
cond3$condition3 <- factor(cond$condition1, levels = c("DOM -> DOM","SUB -> SUB","SUB -> DOM","DOM -> SUB"))
diff.lm4<- glmer(mean_con_ng_ul~period+condition3+ (1|pre_idbatchcage)+ (1|post_idbatch)+(1|batch)+(1|plate), data =cond3,family = Gamma(link = "log"))
summary(diff.lm4)

plot(diff.lm4)
acf(resid(diff.lm4))
qqPlot(resid(diff.lm4))



# Figure 5 -- looking at different domgroup across new cages in just post data
df <- read_csv("manuscript/cort/FullCort.csv") 
head(df)

df<- df %>% dplyr::select(-Well, -final_conc_ng_ul) %>%unique(.)


#changing times
df$time <- ifelse(df$time == "1 hr", "70 min", df$time)
df$time <- ifelse(df$time == "25hr", "25 hr", df$time)
table(df$time)

#domgroup for post data
df$domgroupPost <- ifelse(df$Postrank == 3|df$Postrank==4, "Subordinate", "Dominant")  


df25 <-df %>% filter(group == "reorganized")
df25$type<- ifelse(df25$Prerank == 1, "alphas", "")
df25$type<- ifelse(df25$Prerank == 3, "gammas", df25$type)
df25$type<- ifelse(df25$Prerank == 4, "deltas", df25$type)

df25$type <-factor(df25$type, levels = c("alphas", "gammas", "deltas"))
df25$time <- factor(df25$time, levels = c("70 min", "25 hr"))
df25$domgroupPost<- factor(df25$domgroupPost, levels = c("Dominant", "Subordinate"))

df25 %>% filter(period == "Post") %>% 
  ggplot(., aes(type,mean_con_ng_ul, color =as.factor(domgroupPost)))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21,
                 alpha = 0.8,
                 jitter.height = 0.02, jitter.width = 0.030, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  labs(title = "",
       x = "New Cage",
       y = "Corticosterone (ng/ml)",
       color = "Post Domgroup") +
  scale_color_manual(values = viridis::viridis(4)) +
  scale_fill_manual(values = viridis::viridis(4))+
  ylim(c(0,600))+
  facet_wrap(~time)



df25$period <- as.factor(df25$period)
df25$Prerank <- factor(df25$Prerank, levels = c("1","3", "4"))
df25$Postrank <- factor(df25$Postrank, levels = c("1", "3","4"))
df25$batch <- as.factor(df25$batch)
df25$plate <- as.factor(df25$plate)
df25$pre_idbatchcage <- as.factor(df25$pre_idbatchcage)
df25$condition<- factor(df25$condition, levels = c("control", "same", "ascenders", "descenders"))
df25$post_idbatch <- as.factor(df25$post_idbatch)
df25$group <- as.factor(df25$group)

blah <-  df25 %>% filter(period == "Post")
df25x <- df25 %>% filter(time == "25 hr" & period == "Post")
df25
df1 <- df25 %>%  filter(time == "70 min" & period == "Post")


df1 %>% group_by(type,domgroupPost) %>% 
  summarise(median = median(mean_con_ng_ul),
            lqr = quantile(mean_con_ng_ul,.25),
            uqr = quantile(mean_con_ng_ul,.75))


df25x %>% group_by(type, domgroupPost) %>% 
  summarise(median = median(mean_con_ng_ul),
            lqr = quantile(mean_con_ng_ul,.25),
            uqr = quantile(mean_con_ng_ul,.75))

# linear models for figure 5 - these are terrible models they are not converging idk why something happens we i have domgroupPost and type together in model
str(df25x)
glm <- glmer(mean_con_ng_ul ~ type+ (1|pre_idbatchcage) + (1|pre_id)+ (1|post_id) +(1|batch) + (1|plate), data = df25x, family = Gamma(link= "inverse"))
summary(glm)

glm.dom <- glmer(mean_con_ng_ul ~ domgroupPost+type+(1|pre_idbatchcage) + (1|pre_id)+ (1|post_id) +(1|batch) + (1|plate), data = df25x, family = Gamma(link= "log"))
summary(glm.dom)
plot(glm.dom)

#this works tho? idk 
glm.70min <- glmer(mean_con_ng_ul ~ type+domgroupPost+ type*domgroupPost+(1|pre_idbatchcage) + (1|pre_id)+ (1|post_id) +(1|batch) + (1|plate), data = df1, family = Gamma(link= "log"))
summary(glm.70min)
#but this doesn't converge?
glm.dom70min <- glmer(mean_con_ng_ul ~ domgroupPost+ (1|pre_idbatchcage) + (1|pre_id)+ (1|post_id) +(1|batch) + (1|plate), data = df1, family = Gamma(link= "log"))
summary(glm.dom70min)
plot(glm.70min)

#maybe this will be the best? but this doesn't looked at time 70 min and 25 hr separatly 
glm <- glmer(mean_con_ng_ul ~ domgroupPost+ type+ (1|pre_idbatchcage) + (1|pre_id)+ (1|post_id) +(1|batch) + (1|plate), data = blah, family = Gamma(link= "log"))
summary(glm)
plot(glm)


## Aggression data
#6 Figure - total aggression encounters between domgroupPost in new cages 
x <- read_csv("manuscript/cort/cortAgg70min.csv")
head(x)
x <- x %>% group_by(post_idbatch,Agiven,ARec) %>% 
  mutate(.,tot.ind = as.numeric(Agiven) + as.numeric(ARec))


x$type <- factor(x$type, levels = c("alphas", "gammas", "deltas"))
x$domgroupPost <- ifelse(x$Postrank == 1, "Dominant", "Subordinate")

ggplot(x, aes(tot.ind, mean_con_ng_ul, color = as.factor(domgroupPost)))+
  geom_point(size = 3, alpha = .6)+
  # geom_smooth(method= "lm",se=F)+
  labs(title = "",
       x = "Total Individual Aggression",
       y = "Corticosterone (ng/ml)",
       color = "Post Domgroup") +
  scale_color_manual(values = viridis::viridis(3)) +
  facet_wrap(~type)+
  theme_test()

ggplot(x, aes(mean_con_ng_ul,tot.ind, color = as.factor(domgroupPost)))+
  geom_point(size = 3, alpha = .6)+
  # geom_smooth(method= "lm",se=F)+
  labs(title = "",
       y = "Total Individual Aggression",
       x = "Corticosterone (ng/ml)",
       color = "Post Domgroup") +
  scale_color_manual(values = viridis::viridis(3)) +
  facet_wrap(~type)+
  theme_test()


#linear models 

hist(x$tot.ind) # not normal 

agg.m2 <- glmer(tot.ind ~ mean_con_ng_ul +domgroupPost+type+(1|post_batchcage)+ (1|post_idbatch)+ (1|plate) +(1|batch), data =x, family=Gamma(link="log"))
summary(agg.m2)

