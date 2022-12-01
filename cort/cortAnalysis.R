# cort analysis 
library(lmerTest)
library(lme4)
library(emmeans)
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



### Starting with the control data:
con <- df %>% 
  filter(group == "control") %>% 
  dplyr::select(pre_idbatchcage, period, mean_con_ng_ul, batch, plate, Prerank,domgroup)

head(con)
tail(con)
str(con)


#Just the control data
# cmix <- lmer(cort~period + Prerank + (1|pre_idbatchcage)+ (1|batch)+(1|plate) data =con) 
## maybe we should be doing domgroup

table(con$batch)

table(con$plate, con$batch)
table(con$plate, con$Prerank)

hist(con$mean_con_ng_ul) # not normal 
qqPlot(con$mean_con_ng_ul)


ggplot(con, aes(period,mean_con_ng_ul, group= batch)) +
  geom_point(aes(color = Prerank)) +
  geom_line(aes(group= pre_idbatchcage, color = Prerank))+
  facet_wrap(~batch)

ggplot(con, aes(period,mean_con_ng_ul, group= plate)) +
  geom_point(aes(color = Prerank)) +
  geom_line(aes(group= pre_idbatchcage, color = Prerank))+
  facet_wrap(~plate)


## going to try glmer with family gamma because data is skewed 
head(con)

cgmix <- glmer(mean_con_ng_ul~period + Prerank + (1|pre_idbatchcage)+ (1|batch)+(1|plate), data =con, family = Gamma(link = "inverse"))
summary(cgmix)
AIC(cgmix)

cgmix.l <- glmer(mean_con_ng_ul~period + Prerank + (1|pre_idbatchcage)+ (1|batch)+(1|plate), data =con, family = Gamma(link = "log"))
summary(cgmix.l)
AIC(cgmix.l)

cgmix1 <- glmer(mean_con_ng_ul~period + Prerank + (1|pre_idbatchcage)+ (1|batch)+(1|plate), data =con,family=inverse.gaussian(link = "log"))
summary(cgmix1)
AIC(cgmix1)

cgmix2null <- glmer(mean_con_ng_ul~ 1 + (1|pre_idbatchcage)+ (1|batch)+(1|plate), data =con, family=inverse.gaussian(link = "log"))
summary(cgmix1)
summary(cgmix2null)

AIC(cgmix) 
AIC(cgmix.l)
AIC(cgmix1)# best AIC
AIC(cgmix2null)

anova(cgmix1,cgmix)
anova(cgmix1, cgmix2null)
anova(cgmix1)

xxy <- emmeans(cgmix1, ~ period + Prerank )
xxy
pairs(xxy)

##Assumptions for glmer 
acf(resid(cgmix1)) 
qqPlot(resid(cgmix1))
hist(resid(cgmix1))

plot(x=fitted(cgmix1), y = resid(cgmix1))

durbinWatsonTest(resid(cgmix1))
shapiro.test(resid(cgmix1))


# I think I want to do domgroup instead 


glm.dom <- glmer(mean_con_ng_ul~period + domgroup + (1|pre_idbatchcage)+ (1|batch)+(1|plate), data =con, family = "Gamma")
summary(glm.dom)
AIC(glm.dom)

glm.dom1 <- glmer(mean_con_ng_ul~period + domgroup + (1|pre_idbatchcage)+ (1|batch)+(1|plate), data =con,family=inverse.gaussian(link = "log"))
summary(glm.dom1)
AIC(glm.dom1) ## overall best AIC


##Assumptions for domgroup  glmer, not the best. 
acf(resid(glm.dom1)) 
qqPlot(resid(glm.dom1))
hist(resid(glm.dom1))

plot(x=fitted(glm.dom1), y = resid(glm.dom1))

durbinWatsonTest(resid(glm.dom1))
shapiro.test(resid(glm.dom1))

######################
# Now reorganized data

head(df)
reg <- df %>% 
  filter(group == "reorganized") %>% 
  dplyr::select(pre_idbatchcage, post_idbatch,period,mean_con_ng_ul, batch, plate, Prerank, Postrank, condition,domgroup,domgroupPost)

head(reg)
str(reg)


## checking out a few things
table(reg$condition, reg$plate)
hist(reg$mean_con_ng_ul)  #not as bad as the control, but still leading towards gamma
qqPlot(reg$mean_con_ng_ul)

ggplot(reg, aes(period,mean_con_ng_ul, group= batch)) +
  geom_point(aes(color = Prerank)) +
  geom_line(aes(group= pre_idbatchcage, color = Prerank))+
  facet_wrap(~batch)

ggplot(con, aes(period,mean_con_ng_ul, group= plate)) +
  geom_point(aes(color = Prerank)) +
  geom_line(aes(group= pre_idbatchcage, color = Prerank))+
  facet_wrap(~plate)

## glmers 

reggmix <- glmer(mean_con_ng_ul~period + Prerank + Postrank+ condition+ (1|pre_idbatchcage)+ (1|post_idbatch)+(1|batch)+(1|plate), data =reg, family = Gamma(link = "log"))
summary(reggmix )
AIC(reggmix) #better 

reg.in <- glmer(mean_con_ng_ul~period + Prerank + Postrank+ condition+ (1|pre_idbatchcage)+ (1|post_idbatch)+(1|batch)+(1|plate), data =reg, family = gaussian(link = "inverse"))
summary(reg.in )
AIC(reg.in )

reggmix1 <- glmer(mean_con_ng_ul~period + Prerank + Postrank+ condition+ (1|pre_idbatchcage)+ (1|post_idbatch)+(1|batch)+(1|plate), data =reg,family =inverse.gaussian(link = "log"))
summary(reggmix1)
AIC(reggmix1)

g.nb<- glmer.nb(mean_con_ng_ul~period + Prerank + Postrank+ condition+ (1|pre_idbatchcage)+ (1|post_idbatch)+(1|batch)+(1|plate), data =reg)
summary(g.nb)
AIC(g.nb)


AIC(reggmix) #better 
AIC(reg.in )
AIC(reggmix1)
AIC(g.nb)


# assumptions 
acf(resid(reggmix)) 
qqPlot(resid(reggmix))
hist(resid(reggmix))

plot(x=fitted(reggmix), y = resid(reggmix))

durbinWatsonTest(resid(reggmix))
shapiro.test(resid(reggmix))

## looking at dom group

dom.g <- glmer(mean_con_ng_ul~period + domgroup+ domgroupPost+ condition+ (1|pre_idbatchcage)+ (1|post_idbatch)+(1|batch)+(1|plate), data =reg, family = Gamma(link = "log"))
summary(dom.g)
AIC(dom.g) 

dom.g1 <- glmer(mean_con_ng_ul~period + domgroup+ domgroupPost+ condition+ (1|pre_idbatchcage)+ (1|post_idbatch)+(1|batch)+(1|plate), data =reg,family =inverse.gaussian(link = "log"))
summary(dom.g1)
AIC(dom.g1)

##Assumptions for domgroup  glmer,TERRIBLE 

acf(resid(dom.g)) 
qqPlot(resid(dom.g))
hist(resid(dom.g))

plot(x=fitted(dom.g), y = resid(dom.g))

durbinWatsonTest(resid(dom.g))
shapiro.test(resid(dom.g))


## comparing controls and reorganized 

table(df$group,df$period)
hist(df$mean_con_ng_ul)
qqPlot(df$mean_con_ng_ul)

## just period 1 data and pre data  
head(df)

grp <- df %>%  filter( period!= 2)
## maybe get rid of condition 

head(grp)
grp$condition <- as.factor(grp$condition)

grp$condition1<- ifelse(grp$condition == "ascenders",  "Ascenders", grp$condition)
grp$condition1 <-ifelse(grp$condition1 == "4", "Descenders", grp$condition1)
grp$condition1 <- ifelse(grp$domgroup == "Dominant" & grp$condition1 =="2", "Dominant", grp$condition1)
grp$condition1 <- ifelse( grp$domgroup== "Subordinate" & grp$condition1 == "2", "Subordinate", grp$condition1)
grp$condition1 <- ifelse(grp$condition1 == "1", "Control", grp$condition1)
grp$condition1 <- factor(grp$condition1, levels = c("Control", "Ascenders", "Descenders","Dominant", "Subordinate"))


grp$condition <- factor(grp$condition, levels = c("control","ascenders", "descenders","same"))

## need to figure this out 
glm1 <- glmer(mean_con_ng_ul~group + period +condition1 + Prerank+ Postrank +(1|pre_idbatchcage)+ (1|post_idbatch)+(1|batch)+(1|plate), data =grp, family= Gamma(link ="log"))
summary(glm1)
AIC(glm1)

glm2 <- glmer(mean_con_ng_ul~group + period +condition1 + Prerank+(1|pre_idbatchcage)+ (1|post_idbatch)+(1|batch)+(1|plate), data =grp, family= Gamma(link ="log"))
summary(glm2)
AIC(glm2)

glm3<- glmer(mean_con_ng_ul~condition1+ period+(1|pre_idbatchcage)+ (1|post_idbatch)+(1|batch)+(1|plate), data =grp, family= Gamma(link ="log"))
summary(glm3)
AIC(glm3)


# just control and one hour data
ggplot(grp, aes(group,mean_con_ng_ul, fill =period))+
  geom_boxplot()


## what about both time points?
glm.df<- glmer(mean_con_ng_ul~condition+period +(1|pre_idbatchcage)+ (1|post_idbatch)+(1|batch)+(1|plate), data =df, family= Gamma(link ="log"))
summary(glm.df)
AIC(glm.df)


# all data 
ggplot(df, aes(group,mean_con_ng_ul, fill =period))+
  geom_boxplot()

#NULL model 
glmNull <- glmer(mean_con_ng_ul~1 +(1|pre_idbatchcage)+ (1|post_idbatch)+(1|batch)+(1|plate), data =grp,family = Gamma(link ="log"))
summary(glmNull)
AIC(glmNull)


anova(glm3,glmNull)

tot <- emmeans(glm3, ~condition1 +period)
tot
pairs(tot)

##Assumptions for glmer TERRIBLE  
acf(resid(glm3))  
qqPlot(resid(glm3))
hist(resid(glm3))

plot(x=fitted(glm1), y = resid(glm1))

durbinWatsonTest(resid(glm1))
shapiro.test(resid(glm1))

