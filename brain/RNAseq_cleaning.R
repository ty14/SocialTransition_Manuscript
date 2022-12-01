# Cleaning raw data from Won a little. 
library(tidyverse)

#only did this stuff once:
# one <- read.delim2("manuscript/brain/raw_data/allcounts_firsthalf.txt")
# head(one)
# write.csv(one,"manuscript/brain/raw_data/allcounts_firsthalf.csv",row.names = F)
# 
# two <-read.delim2("manuscript/brain/raw_data/allcounts_secondhalf.txt")
# head(two)
# 
# write.csv(two,"manuscript/brain/raw_data/allcounts_secondhalf.csv", row.names = F)


## joining all the data together 

df <- read_csv("manuscript/brain/raw_data/allcounts_firsthalf.csv")
head(df)
df1 <- read_csv("manuscript/brain/raw_data/allcounts_secondhalf.csv")
head(df1)


gene <- df %>% 
  full_join(df1)
colnames(gene)
write.csv(gene,"manuscript/brain/allcounts.csv", row.names = F)


head(gene)

gene <- read_csv("manuscript/brain/allcounts_fixed.csv")
head(gene)
mpfc <- gene %>%
  select(X, contains("PFC"))

write.csv(mpfc,"manuscript/brain/PFC_counts.csv", row.names = F)

amy<- gene %>%
  select(X, contains("AMY"))

write.csv(amy,"manuscript/brain/AMY_counts.csv", row.names = F)

# think second half of counts had all data in it already? Maybe Won combined everything b4 sending it to me?

# Getting sample information for DESeq2 analysis 

dx <- read_csv("manuscript/brain/raw_data/Sample_info.csv")
head(dx)

# function 
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#id
dx$sample <- gsub(".trim.sam.counts", "", dx$SampleName)
head(dx)

# I sent the core the wrong fucking name B4.PFCB8.2.3 should be B4.PFCB8.2.4
dx$sample <- gsub("B4.PFCB8.2.3", "B4.PFCB8.2.4",dx$sample)
dx$SampleName <- gsub("mPFC", "PFC", dx$SampleName)
dx$SampleName <- gsub("B4.PFCB8.2.3", "B4.PFCB8.2.4",dx$SampleName)
#brain region 
dx$region <- substr(dx$sample,4,6)

#batch 
dx$batch <- substr(dx$sample,8,10)
dx$batch <- sub("\\.\\d+$", "", dx$batch)
dx$batch <- sub("\\.", "", dx$batch)
dx$batch <- sub("B", "", dx$batch)
dx$batch <- paste0("Batch",dx$batch)


##post_id 
dx$post_id <- substrRight(dx$sample,3)
dx$post_id <- gsub("\\.", "-", dx$post_id)

#post_idbatch 

dx$post_idbatch <- paste(dx$post_id,dx$batch)
dx$post_idbatch <- gsub(" ", "", dx$post_idbatch)

head(dx)

dx$batch <- as.numeric(gsub("Batch", "", dx$batch)) 
  
dx <- dx %>%select( -sample, -post_id)

dx1 <- dx

cort <- read_csv("manuscript/cort/FullCort.csv")
head(cort)

cort <- cort %>% dplyr::select(5:7,11:20)
cort


grp <- read_csv("data_raw/groups.csv")
head(grp)

grp$batch <- paste0("Batch", grp$batch)
grp$precage <- paste0("Cage", grp$precage)
grp$pre_idbatch <- paste(grp$preid, grp$batch)
grp$pre_idbatchcage <- paste(grp$pre_idbatch, grp$precage)
grp$pre_idbatchcage <- gsub(" ", "", grp$pre_idbatchcage)

head(grp)
grp <- grp %>% select(4:6,11)

data <- cort %>%  full_join(grp) %>% unique(.) %>% na.omit(.) %>% filter(period =="Post")

## finding aggression given and received 
#libraries 
library(tidyverse)
library(magrittr)
library(lubridate)

##source 
postwl <- read.csv("data_clean/WinLoss_Post.csv", stringsAsFactors = F)
head(postwl)

#assume 25 minutes before lights off was the start time of observations
postwl$obstime <- postwl$ztime + (25*60)
#fix time
postwl$time <- gsub(" hour", 'min', postwl$time)
postwl$time <- gsub("1", '70', postwl$time)
postwl$batch <-paste0("Batch", postwl$batch)

dt <- postwl %>% select(batch, winner, loser, time, obstime)
dt$time
#get observation for the first 70 min of 25 hr reorg
dt <- dt %>% filter(time == " 70min") %>% filter(obstime <4200)
head(dt)
tail(dt)

#split winner and losers 
w <- dt %>% select(batch, winner)
w$post_idbatch <- paste(w$winner, w$batch)
w$post_idbatch <- gsub("Cage ", "", w$post_idbatch)
w$post_idbatch <- gsub(" ", "", w$post_idbatch)
w <- w%>% 
  group_by(winner) %>% 
  mutate(AggGiven70min = n()) %>% 
  unique(.) %>% ungroup(.)
w <- w %>% select(post_idbatch, AggGiven70min)

head(w)
#Aggression Received Post Reorg 70 min
l <- dt %>% select(batch, loser)
l$post_idbatch <- paste(l$loser, l$batch)
l$post_idbatch <- gsub("Cage ", "", l$post_idbatch)
l$post_idbatch <- gsub(" ", "", l$post_idbatch)

l <- l%>% 
  group_by(loser) %>% 
  mutate(AggRec70min = n()) %>% 
  unique(.) %>% ungroup(.) %>% 
  select(post_idbatch,AggRec70min)

head(l)

agg <- w %>% full_join(l) %>%
  mutate_if(is.numeric,coalesce,0)

data70 <- dx1 %>% full_join(agg) %>% na.omit(.) %>% full_join(data)

ds <- read_csv("data_clean/DS_prepost.csv")
head(ds)

ds$post_idbatch <- gsub("Mouse", " ", ds$post_idbatchcage)
ds$post_idbatch <- gsub("Cage", " ", ds$post_idbatch)
ds$post_idbatch <- substr(ds$post_idbatch, 1,11)
ds$post_idbatch <- gsub(" ", "", ds$post_idbatch)
ds$pre_idbatch <- substring(ds$pre_idbatchcage,1,7)
ds$time <- ds$period
x <- ds %>%  select(post_idbatch, pre_idbatch,group,time,Postds,Postrank, Preds, Prerank)

library(data.table)
data70x <-setDT(data70)[setDT(x), on = "post_idbatch", `:=` (pre_idbatch = i.pre_idbatch, group= i.group,time =i.time, Postds =i.Postds, Postrank = i.Postrank,Preds =i.Preds, 
                                                             Prerank = i.Prerank)]

head(data70x)

data70xx <- data70x %>% filter(time!="25hr")

write.csv(data70xx,"manuscript/brain/sample70min_table.csv",row.names = F)
################################################
#25 hr count data and sample_table 

# function 
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

df <- read.delim2("manuscript/brain/raw_data/allcounts_25hr.txt")
head(df)
write.csv(df,"manuscript/brain/raw_data/all_25hcounts.csv", row.names = F)

p <- df %>% 
  dplyr::select(X, contains("P"))
head(p)

write.csv(p,"manuscript/brain/PFC_25hcounts.csv", row.names = F)

a<- df %>%
  select(X, contains("A.b"))
head(a)
#fixing sample mistake
colnames(a)[8] <- "B2.A.b1.2.2.trim.sam.counts"


write.csv(a,"manuscript/brain/AMY_25hcounts.csv", row.names = F)

# getting sample information for 25 hr group.
dx <- read_csv("manuscript/brain/raw_data/Sample_info25hr.csv")
head(dx)

# post_idbatch 2-3Batch1 is suppose to be 2-2Batch1
dx$SampleNames <- gsub('B2.A.b1.2.3.trim.sam.counts', 'B2.A.b1.2.2.trim.sam.counts', dx$SampleNames)

#id
dx$sample <- gsub(".trim.sam.counts", "", dx$SampleNames)
head(dx)

#brain region 
dx$region <- substr(dx$sample,4,5)
dx$region <-gsub("\\.", "",dx$region )

#batch 
dx$batch <- substr(dx$sample,7,8)
dx$batch <- sub("\\.", "", dx$batch)
dx$batch <- paste0("Batch",dx$batch)


##post_id 

dx$post_id <- substrRight(dx$sample,3)
dx$post_id <- gsub("\\.", "-", dx$post_id)

#post_idbatch 

dx$post_idbatch <- paste(dx$post_id,dx$batch)
dx$post_idbatch <- gsub(" ", "", dx$post_idbatch)

dx$batch <- as.numeric(gsub("Batch", "", dx$batch)) 

head(dx)

dx <- dx %>%select( -sample, -post_id)



##Cort and Rank data 

cort <- read_csv("manuscript/cort/FullCort.csv")
head(cort)

cort <- cort %>% dplyr::select(5:7,11:20)
cort


grp <- read_csv("data_raw/groups.csv")
head(grp)

grp$batch <- paste0("Batch", grp$batch)
grp$precage <- paste0("Cage", grp$precage)
grp$pre_idbatch <- paste(grp$preid, grp$batch)
grp$pre_idbatchcage <- paste(grp$pre_idbatch, grp$precage)
grp$pre_idbatchcage <- gsub(" ", "", grp$pre_idbatchcage)

head(grp)
grp <- grp %>% select(4:6,11)


data <- grp %>% full_join(cort)



data <-data %>%  unique(.) %>% filter(period == 'Post') %>% filter(Prerank != 3) %>% filter(Postrank != 3) %>% filter(condition != "control")

#finding missing data # missing data is 3-3Batch1  
head(dx,9)
data %>% filter(post_idbatch == '3-3Batch1')
which(dx$post_idbatch %in% data$post_idbatch)

data <- data[which(data$post_idbatch %in% dx$post_idbatch),] 


## finding aggression given and received 
#libraries 
library(tidyverse)
library(magrittr)
library(lubridate)

##source 
postwl <- read.csv("data_clean/WinLoss_Post.csv", stringsAsFactors = F)
head(postwl)

#assume 25 minutes before lights off was the start time of observations
postwl$obstime <- postwl$ztime + (25*60)
#fix time
postwl$time <- gsub(" hours", 'hr', postwl$time)
postwl$batch <-paste0("Batch", postwl$batch)

dt <- postwl %>% select(batch, winner, loser, time, obstime)

#get observation for the first 70 min of 25 hr reorg
dt <- dt %>% filter(time == "25hr") %>% filter(obstime <4200)
head(dt)
tail(dt)

#split winner and losers 
w <- dt %>% select(batch, winner)
w$post_idbatch <- paste(w$winner, w$batch)
w$post_idbatch <- gsub("Cage ", "", w$post_idbatch)
w$post_idbatch <- gsub(" ", "", w$post_idbatch)
w <- w%>% 
  group_by(winner) %>% 
  mutate(AggGiven70min = n()) %>% 
  unique(.) %>% ungroup(.)
w <- w %>% select(post_idbatch, AggGiven70min)

head(w)
#Aggression Received Post Reorg 70 min
l <- dt %>% select(batch, loser)
l$post_idbatch <- paste(l$loser, l$batch)
l$post_idbatch <- gsub("Cage ", "", l$post_idbatch)
l$post_idbatch <- gsub(" ", "", l$post_idbatch)

l <- l%>% 
  group_by(loser) %>% 
  mutate(AggRec70min = n()) %>% 
  unique(.) %>% ungroup(.) %>% 
  select(post_idbatch,AggRec70min)

head(l)

agg <- w %>% full_join(l) %>%
  mutate_if(is.numeric,coalesce,0)

data25 <- dx %>% full_join(agg) %>% na.omit(.)
data25 <- data25 %>% full_join(data)

ds <- read_csv("data_clean/DS_prepost.csv")
head(ds)

ds$post_idbatch <- gsub("Mouse", " ", ds$post_idbatchcage)
ds$post_idbatch <- gsub("Cage", " ", ds$post_idbatch)
ds$post_idbatch <- substr(ds$post_idbatch, 1,11)
ds$post_idbatch <- gsub(" ", "", ds$post_idbatch)

ds$pre_idbatch <- substring(ds$pre_idbatchcage,1,7)
ds$time <- ds$period
ds <- ds %>%  select(post_idbatch, pre_idbatch,group,time,Postds,Postrank, Preds, Prerank)

library(data.table)
data25x <-setDT(data25)[setDT(ds), on = "post_idbatch", `:=` (pre_idbatch = i.pre_idbatch, group= i.group,time =i.time, Postds =i.Postds, Postrank = i.Postrank,Preds =i.Preds, 
                                                             Prerank = i.Prerank)]

head(data25x)


write.csv(data25x,"manuscript/brain/sample25hr_table.csv",row.names = F)
