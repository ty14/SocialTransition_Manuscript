### Compare Day 1 to Day 2 

library(tidyverse)

xf <- read_csv("manuscript/behavior/Post_WL.csv")
xf <- xf[-407,]   #batch 10, typo of ids.

head(xf)
table(xf$time)

xf %>% 
  filter(time == "25 hours") %>%
  filter(group == "reorganized") %>%
  group_by(post_batchcage) %>% 
  mutate(day = as.numeric(date - min(date) + 1)) -> xfd

head(xfd)

# create list of dataframes for each individual cage
l <- split(xfd, list(xfd$post_batchcage, xfd$day))


## check if alphas are the same.

# just keep winner loser columns
l.wl <- l %>% map(~ ungroup(.)) %>% map(~ select(., winner,loser) %>% as.data.frame)

#check DS
l.ds <- lapply(l.wl, function(x) compete::ds(compete::get_wl_matrix(x)))

# all alphas are the same on both days
matrix(unlist(lapply(l.ds, function(x) names(rev(x[order(x)])[1]))),ncol = 2)


## or, did individuals of rank 1 lose to any of rank 2/3/4 on day 2?
# did individuals of rank 2 lose to any of rank 3/4 on day 2?

#Get ranks day 1... and use same matrix to enter wins/losses

l.mat <- lapply(l.wl, function(x) reshape2::melt(compete::org_matrix(compete::get_wl_matrix(x), method="ds")))
l.mat0 <- lapply(l.wl, function(x) reshape2::melt(compete::org_matrix(compete::get_di_matrix(compete::get_wl_matrix(x)), method="ds")))


substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

### Raw matrix
Map(cbind, 
    l.mat, 
    ids = names(l.mat), 
    day = substrRight(names(l.mat),1), 
    cage = substr(substrRight(names(l.mat),3),1,1),
    batch = substr(gsub("Batch ", "", names(l.mat)),1,1)
) -> l.mat1

df1 <- data.table::rbindlist(l.mat1[1:24])
df2 <- data.table::rbindlist(l.mat1[25:48])
df1$ids<-df2$ids<-NULL
colnames(df2)[3]<-"value2"
dfx <- full_join(df1 %>% select(-day),df2 %>% select(-day))

dfx


### Binary matrix
Map(cbind, 
    l.mat0, 
    ids = names(l.mat0), 
    day = substrRight(names(l.mat0),1), 
    cage = substr(substrRight(names(l.mat0),3),1,1),
    batch = substr(gsub("Batch ", "", names(l.mat0)),1,1)
) -> l.mat0

df01 <- data.table::rbindlist(l.mat0[1:24])
df02 <- data.table::rbindlist(l.mat0[25:48])
df01$ids<-df02$ids<-NULL
colnames(df02)[3]<-"value2"
df0x <- full_join(df01 %>% select(-day),df02 %>% select(-day))

df0x$value2 <- ifelse(is.na(df0x$value2),0,df0x$value2)
df0x

# add in ranks.
df0x$wrank <- rep(1:4,4)
df0x$lrank <- rep(1:4,each=4)
df0x

# Alphas
df0x %>% filter(wrank==1, lrank>1)
df0x %>% filter(wrank>1, lrank==1)

#Betas
df0xb <- df0x %>% filter(wrank>2, lrank==2)

#Batch 7 cage 3, Cage 3-1 Cage 5-1  [1 fight]
#Batch 4 cage 4,  Cage 4-3 Cage 3-4 [1 fight]
#Batch 2 cage 2,  Cage 5-4 Cage 2-3

table(df0xb$cage,df0xb$batch)
binom.test(21,24)


#Gammas
df0x %>% filter(wrank>3, lrank==3)  #cage 3 batch 4;  Cage 3-1 .. Cage 4-4 [1 fight]

l.mat
