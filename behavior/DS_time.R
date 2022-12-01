library(tidyverse)

#source

#####Pre data

df <- read_csv("manuscript/behavior/Pre_WL.csv")
head(df)

# Get rid of ties - David's Scores won't allow
df <- df[help$score!=0.5,]


### Just do on one cage


x <- df[df$pre_batchcage=="Batch 8 Cage 2",]
table(x$winner,x$loser)
head(x)
x <- x %>% arrange(ztime)
x

table(x$score)
table(df$score)

results <- NULL
for(i in 1:nrow(x)){
  results[[i]] <- ds(get_wl_matrix(x[1:i,2:3])) # rows 1 thru i
}


ids <- unlist(lapply(results,names))
vals <- unlist(results)
dx <- data.frame(ids,vals, time = rep(1:nrow(x), times =     unlist(lapply(results, length)) ) )
head(dx)     

dx[1:20,]  

# just to look
ggplot(dx, aes(x=time, y=vals, color=ids)) + geom_line(lwd=1) + theme_minimal()

#


### Convert to a function.

ds_df <- function(x){
  
  x <- x %>% arrange(ztime)
  
  results <- NULL
  for(i in 1:nrow(x)){
    results[[i]] <- ds(get_wl_matrix(x[1:i,2:3])) # rows 1 thru i
  }
  
  ids <- unlist(lapply(results,names))
  vals <- unlist(results)
  dx <- data.frame(ids,vals, time = rep(1:nrow(x), times =     unlist(lapply(results, length)) ) )
  
  return(dx)
}


## 
head(help)

help.sp <- split(help, help$pre_batchcage)

help.sp.dsdf <- lapply(help.sp, ds_df)

help.sp.dsdf[[4]]

######

ggplot(help.sp.dsdf[[54]], aes(x=time, y=vals, color=ids)) + geom_line(lwd=1) + theme_minimal()
ggplot(help.sp.dsdf[[22]], aes(x=time, y=vals, color=ids)) + geom_line(lwd=1) + theme_minimal()


######

## One chart summarizing all the David's Scores.


names(help.sp.dsdf)

dy <- do.call('rbind', Map(cbind, help.sp.dsdf, pre_batchcage = names(help.sp.dsdf)))
dy$unique_id <- paste(dy$pre_batchcage, dy$ids, sep="-")

head(dy)


dy.summary <- dy %>% group_by(pre_batchcage) %>% filter(time==max(time)) %>% mutate(rank = dense_rank(-vals))

dy$rank <-  dy.summary$rank[match(dy$unique_id,dy.summary$unique_id)]

head(dy)
tail(dy)

dy.ds <-dy %>% ungroup() %>%
  group_by(rank,time) %>%
  summarize(meanx = mean(vals),
            sdx = sd(vals),
            n= n(),
            medianx = median(vals)) %>%
  mutate(sem = sdx/sqrt(n)) %>% 
  filter(!is.na(sem)) %>% 
  mutate( CI_lower_mean = meanx + qt((1-0.95)/2, n - 1) * sem,
          CI_upper_mean = meanx - qt((1-0.95)/2, n - 1) * sem) %>%
  mutate( CI_lower_med = medianx + qt((1-0.95)/2, n - 1) * sem,
          CI_upper_med = medianx - qt((1-0.95)/2, n - 1) * sem) %>%
  ungroup()

head(dy.ds)

str(dy.ds)
tail(dy.ds)

dy.ds$rank <- as.factor(dy.ds$rank)

library(viridis)

ds_mean_pre <- ggplot(dy.ds, aes(x=time,y=meanx,color=rank)) + 
  geom_line(lwd=1) +
  geom_ribbon(aes(ymin=CI_lower_mean,ymax=CI_upper_mean,fill=rank),alpha=0.4) +
  ylab("Mean David Score") +
  xlab("# of Aggressive Interactions")+
  ggtitle("Mean David Score during Hierarchy Formation")+
  theme_minimal()

ds_mean_pre155 <- dy.ds %>% 
  filter(time < 155) %>% 
  ggplot(., aes(x=time,y=meanx,color=rank)) + 
  geom_line(lwd=1) +
  geom_ribbon(aes(ymin=CI_lower_mean,ymax=CI_upper_mean,fill=rank),alpha=0.4) +
  ylab("Mean David Score") +
  xlab("# of Aggressive Interactions")+
  # ggtitle("Mean David Score during Hierarchy Formation")+
  theme_minimal()


ds_mean_pre <- ggplot(dy.ds, aes(x=time,y=medianx,color=factor(rank))) + 
  geom_line(lwd=1) + 
  geom_ribbon(aes(ymin=CI_lower_med, ymax=CI_upper_med,fill=factor(rank)),alpha=0.4) +
  ggtitle("Median David Score of Ranks during Hierarchy Formation")+
  ylab("Median David Score") +
  xlab("# of Aggressive Interactions")+
  theme_minimal()

ds_med_pre155 <- dy.ds %>% 
  filter(time < 155) %>% 
  ggplot(., aes(x=time,y=medianx,color=factor(rank))) + 
  geom_line(lwd=1) + 
  geom_ribbon(aes(ymin=CI_lower_med, ymax=CI_upper_med,fill=factor(rank)),alpha=0.4) +
  ggtitle("Median David Score of Ranks during Hierarchy Formation")+
  ylab("Median David Score") +
  xlab("# of Aggressive Interactions")+
  theme_minimal()