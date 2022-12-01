#CORT DATA
## missing data from Batch 1 - mice 4-1 (control), 4-3(control), 3-3(sameReorg 4-4)
### * pre and post missing, looking back I probably should have at least done the post. 
library(tidyverse)


df <- read_csv("data_raw/CORTraw.csv")
head(df)
tail(df)

# function to get rid of NA in Well ID column 
fill_well_id <- function(df){
  for(i in 1:nrow(df)){
    ifelse(is.na(df$Well_ID[i]),
           df$Well_ID[i]<-df$Well_ID[i-1],
           df$Well_ID[i])
  }
  return(df)
}


df1 <- fill_well_id(df)

colnames(df1)[5] <- "OD"



head(df1)

# calculate %B/B0 for all reading first
get_B_B0_perc <- function(df){
  df %>% 
    filter(grepl('NSB', Well_ID)) %>% 
    summarize(B0 = mean(OD)) %>% 
    unlist()-> NSB
  
  df1 <- df %>% 
    mutate(net_OD = OD - NSB)
  
  df1 %>% 
    filter(grepl('BLK', Well_ID)) %>% 
    summarize(B0 = mean(net_OD)) %>% 
    unlist()-> B0
  
  df2 <- df1 %>% 
    mutate(B_B0_perc = net_OD/B0*100) %>% 
    filter(!grepl('NSB', Well_ID)) %>% 
    filter(!grepl('BLK', Well_ID)) 
  
  return(df2)
}


get_B_B0_perc(df1)

# data frame listing concentration & well ID; same throughout batches

cort_list <- split(df1, df1$plate)
cort_B_B0_perc_list <- cort_list %>% map(get_B_B0_perc)

conc_df = data.frame(
  
  Standard  = cort_list[[1]] %>% 
    filter(grepl('STD',Well_ID)) %>% 
    select(Well_ID) %>% 
    unique(),
  
  conc = 5000/(2^(0:8))
)

conc_df


# Standard Curve for each data frame (each plate)

# According to the manual I have to fit 4PL (4-parameter logistic function)
# stole code from: https://stats.stackexchange.com/questions/61144/how-to-do-4-parametric-regression-for-elisa-data-in-r


fit_4PL_standard_curve <- function(df){
  library('drc')
  std_df <- df %>% 
    filter(grepl('STD',Well_ID)) %>% 
    left_join(conc_df)
  
  std_mod <- drm(B_B0_perc ~ conc, data = std_df , 
                 fct = LL.4(names=c("Slope","Lower Limit","Upper Limit","ED50")))
  
  return(std_mod)
}

std_model_list <- cort_B_B0_perc_list %>% map(fit_4PL_standard_curve) 


dev.off()
for (i in 1:length(std_model_list)){
  plot(std_model_list[[i]], main = names(std_model_list[i]))
}


for (i in 1:length(std_model_list)){
  print(names(std_model_list)[i])
  print(summary(std_model_list[[i]]))
  
}



# once done inspecting model and confirmed each model is correct,


get_X_4PL_curve <- function(model, Y) {
  # Y = (A-D)/(1+(X/C)^B) + D, where Y = %B/B0, X = concentration
  
  A = model$coefficients[3]# A = Upper limit
  D = model$coefficients[2]# D = Lower Limit
  B = model$coefficients[1] # B = Slope
  C = model$coefficients[4] # C = ED50
  
  # So, to get X based on Y, 
  
  X = C*((A-D)/(Y-D)-1)^(1/B) 
  
  return(X)
}

# I 'could have' written code wihtout forloop but oh well, :P
cort_conc_list <- list()

for(i in 1:length(cort_B_B0_perc_list)){
  
  cort_conc_list[[i]] <- cort_B_B0_perc_list[[i]] %>% 
    filter(!grepl('STD',Well_ID)) %>% 
    mutate(conc = get_X_4PL_curve(std_model_list[[i]],B_B0_perc))
  names(cort_conc_list)[i] = names(cort_B_B0_perc_list[i])
}


cort_conc_list %>% map(head)

conc_df <- cort_conc_list %>% map2_dfr(.,names(.), ~mutate(.x, batchID = .y))
head(conc_df)
tail(conc_df)

library(tidyverse)
# Adding dilution factor!
final_cort_conc <- conc_df %>% 
  mutate(dilution_factor = 100)%>% 
  mutate(final_conc_ng_ul = conc * dilution_factor/1000)

df1<- final_cort_conc %>% dplyr::select(Name,Well,plate,dilution_factor,final_conc_ng_ul)

## getting rid of NAs and making subject names for every row
library(zoo)
df1 <-na.locf(df1)

head(df1)

## getting mean final concentration

avg<- df1 %>% group_by(Name) %>% summarise(mean_con_ng_ul = mean(final_conc_ng_ul)) 

df1 <- df1 %>% full_join(avg)

head(df1)


## creating columns for idenifying information and to join with David score data 

df1$batch <- as.numeric(str_sub(df1$Name, 2,3))
df1$batch <- paste0("Batch", df1$batch)
head(df1)
tail(df1)

df1$period <- ifelse(grepl("pre", df1$Name), "Pre", "Post")
head(df1)

df1$post_id <- str_sub(df1$Name,4,8)
df1$post_id  <- gsub(' p',"", df1$post_id )
df1$post_id <- gsub(" ","", df1$post_id )
head(df1)
tail(df1)


df1$precage <- str_sub(df1$post_id , 1,1)
df1$precage<- paste0("Cage", df1$precage)
df1$pre_id <- str_sub(df1$post_id,3,3)

df1$post_idbatch <- paste(df1$post_id, df1$batch)
df1$post_idbatch<- gsub(" ","", df1$post_idbatch)

df1$pre_idbatch <- paste(df1$pre_id, df1$batch)

df1$pre_idbatchcage<- paste(df1$pre_idbatch, df1$precage)
df1$pre_idbatchcage <- gsub(" ","", df1$pre_idbatchcage)

head(df1)

df1$batch <- gsub("\\D+", "", df1$batch)


head(df1)
tail(df1)

df1 <- df1 %>% dplyr::select(-Name)

head(df1)

write.csv(df1,"manuscript/cort/CleanCort.csv",row.names = F)

####### Now adding in the ds and condition data 

ds <- read_csv("data_clean/DS_prepost.csv")
head(ds)


ds$post_idbatch <-str_sub(ds$post_idbatchcage, 6,15)
ds$post_idbatch <- gsub("C", "", ds$post_idbatch)

ds$pre_idbatch <-str_sub(ds$pre_idbatchcage, 1,8)
ds$pre_idbatch <- gsub("C", "", ds$pre_idbatch)


colnames(ds)

ds <- ds %>% dplyr::select(batch,group,period,Postds, Postrank, Preds,Prerank,Precage, PreID, post_idbatch,pre_idbatch, pre_idbatchcage)
head(ds)

colnames(ds)[3] <- "time"
colnames(ds)[8] <- "precage"
colnames(ds)[9] <- "pre_id"

ds$precage <- gsub(" ", "", ds$precage)
df1$pre_idbatch <- gsub(" ", "", df1$pre_idbatch)
# now joining data together
head(ds)
head(df1)

colnames(ds)
colnames(df1)

str(ds)
str(df1)
df1$batch <- as.numeric(df1$batch)
ds$pre_id <- as.character(ds$pre_id)


xx <- df1 %>% full_join(ds)
## probably should not need to do this, but idk why it works and gets me back to 502. 
## where are the 139 nas coming from the rank 2s?
dfx <- na.omit(xx)

## sanity check 
unique(dfx$pre_idbatchcage, df1$pre_idbatchcage)
dfx$pre_idbatchcage %in% df1$pre_idbatchcage

## looks good, but where are the nas in xx coming from? -_-


## now need to add in all the conditions
con <- read_csv("data_clean/CortElisa1.csv")
con1 <- read_csv("data_clean/CortElisa2.csv")


el <- con %>% rbind(con1)

el$post_idbatch <- paste(el$id, el$batch)
head(el)

el$post_idbatch <- gsub(" ", "", el$post_idbatch)
head(el)

el$batch <- as.numeric(gsub("Batch", "", el$batch))

head(el)
colnames(el)[2] <- "post_id"
colnames(el)[7] <- "time"

full_cort <- dfx %>% 
  full_join(el)

## filter out the the 3 we don't have because I couldn't find them 

dc <- na.omit(full_cort) ## great at 502 which is correct

### need to change the controls that are 2 and 3 to 1 and 4

table(dc$Postrank)
table(dc$Postrank, dc$condition)
str(dc)

dc$Postrank <- ifelse(dc$condition== "control" & dc$Postrank == 2 & dc$Prerank == 4, 4, dc$Postrank)
dc$Postrank <- ifelse(dc$condition== "control" & dc$Postrank == 2 & dc$Prerank == 3, 3, dc$Postrank)

dc$Postrank <- ifelse(dc$condition== "control" & dc$Postrank == 3 & dc$Prerank == 4, 4, dc$Postrank)
dc$Postrank <- ifelse(dc$condition== "control" & dc$Postrank == 4 & dc$Prerank == 3, 3, dc$Postrank)

check <- dc %>% filter(condition == "control")

table(dc$Postrank)
table(dc$Postrank, dc$condition) #seems correct now
table(dc$Prerank, dc$Postrank)

head(dc)



write.csv(dc,"manuscript/cort/FullCort.csv",row.names = F)



dc %>% 
  ggplot(aes(x = mean_con_ng_ul))+
  geom_histogram()+
  facet_wrap(~plate)+
  theme_bw()

#yikes might be batch effects. Plate 5 is all low -_-

dc %>% 
  ggplot(aes(as.factor(plate), mean_con_ng_ul))+
  geom_boxplot(alpha = 0.2, outlier.color = NA)+
  geom_jitter(width = 0.1)+
  theme_bw() 



########################
#graphs 

### first graph boxplot of CORT (yaxis)  vs  group (to-be-reorganized  vs control on x-axis)
library(viridis)

g1 <- dc %>%filter(period == "Pre")
  
g1$group <- ifelse(g1$group == "reorganized", "To-Be-Reorganized","Control")

  
  ggplot(g1, aes(group, mean_con_ng_ul,color =group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21,
                 alpha = 0.8,
                 jitter.height = 0.02, jitter.width = 0.030, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  labs(title = "Corticosterone Proir to Reorganization",
    x = "",
       y = "Corticosterone (ng/ml)") +
    scale_color_manual(values = viridis::viridis(4)) +
    scale_fill_manual(values = viridis::viridis(4))+
    newggtheme+
    theme( legend.position = 'none')
  
  g1 %>% group_by(group) %>%
    summarise(median = median(mean_con_ng_ul),
              lqr = quantile(mean_con_ng_ul,.25),
              uqr = quantile(mean_con_ng_ul,.75))
  
    


### second graph boxplot of CORT (yaxis)  vs rank (1,2,3,4  regardless of whether from to-be-reorganized group or control group).
g1$Prerank <- as.factor(g1$Prerank)
    
  
  ggplot(g1,aes(Prerank,mean_con_ng_ul,color, color=Prerank))+
    geom_boxjitter(outlier.color = NA, jitter.shape = 21,
                   alpha = 0.8,
                   jitter.height = 0.02, jitter.width = 0.030, errorbar.draw = TRUE,
                   position = position_dodge(0.85)) +
      labs(title = "Corticosterone Proir to Reorganization",
        x = "Rank",
           y = "Corticosterone (ng/ml)") +
  scale_color_manual(values = viridis::viridis(4)) +
    scale_fill_manual(values = viridis::viridis(4))+
    newggtheme+
    theme( legend.position = 'none')

  
g1$domgroup <- ifelse(g1$Prerank == 3|g1$Prerank==4, "Subordinate", "Dominant")  
table(g1$Prerank,g1$domgroup) # seems werid to have only 74 rank 3 and rank 4, but 100 rank 1 


## second graph with domgroups instead
dompre <-ggplot(g1,aes(domgroup,mean_con_ng_ul,color, color=domgroup))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21,
                 alpha = 0.8,
                 jitter.height = 0.02, jitter.width = 0.030, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  labs(title = " Proir to Reorganization",
    x = "",
       y = "Corticosterone (ng/ml)") +
  scale_color_manual(values = viridis::viridis(4)) +
  scale_fill_manual(values = viridis::viridis(4))+
  newggtheme+
  theme( legend.position = 'none')


### reorganized graphs.

# boxplot of CORT (yaxis)  vs  group (to-be-reorganized  vs control on x-axis), add separate boxplots based on time taken  (pre, 1hr, 25hr)

dc$period <- factor(dc$period, levels = c("Pre", "Post"))
dc$group <- ifelse(dc$group == "control", "Control", "Reorganized")


dc$Prerank <- as.factor(dc$Prerank)

dc$domgroup_pre <- ifelse(dc$Prerank == 3|dc$Prerank==4, "Subordinate", "Dominant")  


dc %>% 
  ggplot(., aes(x=period , y = mean_con_ng_ul)) +
  geom_point(aes(color =Prerank), alpha = .6, size = 3) +
  geom_line(aes(group = pre_idbatchcage, color =Prerank), size = .5) +
  ylim(0,500)+
  labs(x = "",
       y = "Corticosterone (ng/ml)",
       legend= "Group") +
    scale_color_manual(name = "Pre Rank", values = viridis::viridis(3)) +
  facet_grid(time~group, space="free", scales="free")+
  newggtheme 

## second graph

dc$domgroup <- ifelse(dc$Postrank == 3|dc$Postrank==4, "Subordinate", "Dominant")  
table(dc$Postrank,dc$domgroup)

ggplot(dc, aes(period, mean_con_ng_ul))+
  geom_point(aes(color =time), alpha = .6, size = 3) +
  geom_line(aes(group = pre_idbatchcage, color =time), size = .5) +
  labs(x = "",
       y = "Corticosterone (ng/ml)") +
  scale_color_manual(values = viridis::viridis(4)) +
  facet_wrap(~domgroup)+
  newggtheme 


## third.
 dompost <- dc %>% filter(period == "Post") %>% 
ggplot(.,aes(domgroup,mean_con_ng_ul,color, color=domgroup))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21,
                 alpha = 0.8,
                 jitter.height = 0.02, jitter.width = 0.030, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  labs(title = "Post Reorganization",
       x = "",
       y = "Corticosterone (ng/ml)") +
  scale_color_manual(values = viridis::viridis(4)) +
  scale_fill_manual(values = viridis::viridis(4))+
  newggtheme+
  theme( legend.position = 'none')

library(gridExtra)

grid.arrange(dompre,dompost, nrow=1)


##4 panel figure of the slope graphs, with the columns being dom/sub and the rows being 1hr/25 hr

## changing time to Control, 1 hr == 70 min, and 25 hour


dc$time <- ifelse(dc$group == "Control", "Control", dc$time)
dc$time <- ifelse(dc$time == "1 hr", "70 min", dc$time)
dc$time <- ifelse(dc$time == "25hr", "25 hr", dc$time)
dc$time <- factor(dc$time, level = c("70 min", "25 hr", "Control"))

g2 <- dc %>% group_by(domgroup_pre,time,period) %>%
  summarize(meanx = mean(mean_con_ng_ul),
            sdx = sd(mean_con_ng_ul),
            n= n(),
            medianx = median(mean_con_ng_ul)) %>%
  mutate(sem = sdx/sqrt(n)) %>% 
  filter(!is.na(sem)) %>% 
  mutate( CI_lower_mean = meanx + qt((1-0.95)/2, n - 1) * sem,
          CI_upper_mean = meanx - qt((1-0.95)/2, n - 1) * sem) %>%
  mutate( CI_lower_med = medianx + qt((1-0.95)/2, n - 1) * sem,
          CI_upper_med = medianx - qt((1-0.95)/2, n - 1) * sem) %>%
  ungroup()

g2 <- dc %>% dplyr::select(pre_idbatchcage,domgroup, domgroup_pre,period, time,mean_con_ng_ul) %>% full_join(g2)

post <- dc %>% 
  group_by(domgroup,time,period) %>%
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


post <-dc %>% dplyr::select(pre_idbatchcage,domgroup, domgroup_pre,period, time,mean_con_ng_ul) %>% full_join(post)



  ggplot(g2, aes(x=period , y = mean_con_ng_ul)) +
    geom_ribbon(aes(ymin = CI_lower_med,
                    ymax = CI_upper_med, group=domgroup_pre,fill=domgroup_pre),alpha=.9) +
    geom_line(aes(y=medianx, group = domgroup_pre, color = domgroup_pre), size=1.5)+
    geom_line(aes(group = pre_idbatchcage), color= "gray", size = .4, alpha=.4) +
    geom_point(aes(group =domgroup_pre), color ="gray", alpha = .4, size = 1)+
  ylim(0, 500)+
  labs(   x = "",
      y=  "Corticosterone (ng/ml)") +
  scale_color_manual(values = viridis::viridis(3)) +
  facet_grid(time~domgroup_pre, space="free", scales="free")+
  newggtheme +
    theme(legend.position = "none")
  
  
dc %>% group_by(domgroup_pre,time,period) %>%
    summarise(median = median(mean_con_ng_ul),
              lqr = quantile(mean_con_ng_ul,.25),
              uqr = quantile(mean_con_ng_ul,.75))


## what about if you group by domgroup_post which probably makes more sense.. 


ggplot(post, aes(x=period , y = mean_con_ng_ul)) +
  geom_ribbon(aes(ymin = lower_medp,
                  ymax = upper_medp, group=domgroup,fill=domgroup),alpha=.9) +
  scale_fill_discrete(guide=FALSE)+
   geom_line(aes(group = pre_idbatchcage),color ="gray", size = .4, alpha=.4) +
   geom_point(aes(group =domgroup),color="gray", alpha = .4, size = 1)+
  geom_line(aes(y=median_post, group = domgroup, color = domgroup), size=1.5)+
  ylim(0, 500)+
  labs(x = "",
      y=  "Corticosterone (ng/ml)",
      color = "Post Domgroup") +
  scale_color_manual(values = viridis::viridis(3)) +
  facet_grid(time~domgroup_pre, space="free", scales="free")+
  newggtheme




ggplot(post, aes(x=period , y = mean_con_ng_ul)) +
  geom_ribbon(aes(ymin = lower_meanp,
                  ymax = upper_meanp, group=domgroup,fill=domgroup),alpha=.9) +
  scale_fill_discrete(guide=FALSE)+
  geom_line(aes(group = pre_idbatchcage),color ="gray", size = .4, alpha=.4) +
  geom_point(aes(group =domgroup),color="gray", alpha = .4, size = 1)+
  geom_line(aes(y=mean_post, group = domgroup, color = domgroup), size=1.5)+
  ylim(0, 500)+
  labs(x = "",
       y=  "Corticosterone (ng/ml)",
       color = "Post Domgroup") +
  scale_color_manual(values = viridis::viridis(3)) +
  facet_grid(time~domgroup_pre, space="free", scales="free")+
  newggtheme



### thinking about 1 hr data for dom group. 

one1 <- dc %>% filter(time !="25 hr") 
one <- dc %>% filter(time == "70 min") %>% dplyr::select(pre_idbatchcage,domgroup, domgroup_pre,period, time,mean_con_ng_ul,condition)
  
  
 sum <-  dc %>% filter(time == "70 min") %>% 
group_by(domgroup_pre,time,period,condition) %>%
  summarize(sdx = sd(mean_con_ng_ul),
            n= n(),
            medianx = median(mean_con_ng_ul)) %>%
  mutate(semx = sdx/sqrt(n)) %>% 
  filter(!is.na(semx)) %>%
  mutate( lower_med = medianx + qt((1-0.95)/2, n - 1) * semx,
          upper_med = medianx - qt((1-0.95)/2, n - 1) * semx) %>%
  ungroup()

one <- one %>%full_join(sum)


ggplot(one, aes(period,mean_con_ng_ul))+
  geom_ribbon(aes(ymin = lower_med,
                  ymax = upper_med, group=domgroup,fill=domgroup),alpha=.9) +
  scale_fill_discrete(guide=FALSE)+
  geom_line(aes(group = pre_idbatchcage),color ="gray", size = .4, alpha=.4) +
  geom_point(aes(group =domgroup),color="gray", alpha = .4, size = 1)+
  geom_line(aes(y=medianx, group = domgroup, color = domgroup), size=1.5)+
  facet_wrap(~condition)+
  labs(x = "",
      y=  "Corticosterone (ng/ml)",
     color = "Post Domgroup") +
  scale_color_manual(values = viridis::viridis(3))+
  newggtheme


