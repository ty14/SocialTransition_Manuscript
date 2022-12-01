## just get clean win /loss data 

library(tidyverse)

## Read in Data

temp <- list.files(path="data_raw/",pattern="*.csv")
xfiles <- lapply(temp, function(x) read.csv(paste0("data_raw/",x)) )                 
lapply(xfiles, head)
lapply(xfiles, colnames)

prepost <- ifelse(grepl("Pre", temp), "Pre", "Post")


temp[prepost=="Pre"]
pre.l <- xfiles[prepost=="Pre"]
lapply(pre.l, head)


xfiles<- pre.l[c(1:14)]

lapply(xfiles, head)

lapply(xfiles, tail)


# Read in behaviors file
behavs <- read.csv("code_carpentry/behavs.csv")
head(behavs)


# Check have all behaviors

## Get all Behavioral Combinations
behavsX <-
  unique(
    unlist(
      lapply(xfiles, 
             function(x) 
               paste(x$Animal.1..Behavior, x$Animal.2..Behavior, sep=" --- "))
    )
  )



behavsX %in% behavs$behavior # if have any FALSE need to update 'behavs' file.



# Add in score column.... X = scores in behavs 
add_score <- function(df){
  df$behav <-  paste(df$Animal.1..Behavior, df$Animal.2..Behavior, sep=" --- ")
  df$score <- behavs$X[match(df$behav, behavs$behavior)]
  return(df)
}



xfiles <- lapply(xfiles, add_score)
lapply(xfiles, head)



# step 1.  Get Timestamp correct.
# note some datasets have seconds, others do not.
get_times <- function(dataf){
  dataf$ts <- strptime(as.character(dataf$Timestamp),'%m/%d/%Y %H:%M:%S')
  dataf <- dataf[order(dataf$ts),]
  return(dataf)
}

xfiles <- lapply(xfiles, get_times)
lapply(xfiles, head)


# Remove Starts,Ends,NA, when id1==id2  # Remove Starts,Ends,NA, when id1==id2
clean_data <- function(df){
  df <- df[!is.na(df$score),]
  df <- df[df$Start.End!="Start",]
  df <- df[df$Start.End!="End",]
  df <- df[df$Animal.1!=df$Animal.2,]
}

xfiles <- lapply(xfiles, clean_data)

lapply(xfiles, head)

# get winner/loser
get_wl <- function(df){
  df$winner <- ifelse(df$score==0, df$Animal.2, df$Animal.1)
  df$loser <- ifelse(df$score==0, df$Animal.1, df$Animal.2)
  return(df)
}

xfiles <- lapply(xfiles, get_wl)
lapply(xfiles, head)


# only keep time, cage, winner, loser, score
lapply(xfiles, colnames)
xfiles1 <- lapply(xfiles, function(x) x[,c('ts','Cage','winner','loser','score')])
lapply(xfiles1, head)
lapply(xfiles1, tail)





### Add Batch
batchx <-  c("1", "10","11", "12", "13", "14", "2", "3", "4","5", "6", "7", "8", "9")


# batch <- substr(temp, 1, 11)
myfiles1 <- Map(cbind, xfiles1, batch = batchx)

lapply(myfiles1, head)
lapply(myfiles1, tail)


# Fix score
# doing in loop as doesn't seem to work with lapply, don't know why
## I am not sure what this does ?
for(i in 1:length(myfiles1)){
  myfiles1[[i]]$score <- ifelse(myfiles1[[i]]$score==0.5, 0.5, 1)
}

head(myfiles1[[1]])


##  somethings I copy code and hope it works. 

## just getting pre data 
WL_df <- do.call('rbind',myfiles1)
WL_df$batch.cage <- paste(WL_df$batch,WL_df$Cage)
WL_df$batch.cage<- paste0("Batch ", WL_df$batch.cage)

head(WL_df)
tail(WL_df)

#remove batches/cages not required
WL_df <- WL_df[WL_df$batch.cage !="Batch 10 Cage 2",]
WL_df <- WL_df[WL_df$batch.cage !="Batch 8 Cage 5",]
WL_df <- WL_df[WL_df$batch.cage !="Batch 7 Cage 2",]
WL_df <- WL_df[WL_df$batch.cage !="Batch 6 Cage 5",]
WL_df <- WL_df[WL_df$batch.cage !="Batch 3 Cage 3",]
WL_df <- WL_df[WL_df$batch.cage !="Batch 3 Cage 5",]

head(WL_df)

unique(WL_df$batch.cage)



## adding in the controls 
controls <- c("Batch 1 Cage 4","Batch 2 Cage 4","Batch 3 Cage 1","Batch 3 Cage 2","Batch 3 Cage 4",
              "Batch 4 Cage 1","Batch 5 Cage 1","Batch 9 Cage 2","Batch 11 Cage 5",
              "Batch 12 Cage 3","Batch 13 Cage 5","Batch 14 Cage 3","Batch 14 Cage 5")


WL_df$group <- ifelse(WL_df$batch.cage %in% controls, "control", "reorganized")  

table(WL_df$group)


## adding in the timepoints

onehr <- c("8", "9","10", "11", "12", "13", "14")

WL_df$time <- ifelse(WL_df$batch %in% onehr, " 1 hour", "25 hours")

WL_df$ts <- as.POSIXct(WL_df$ts)

head(WL_df)


help <- WL_df %>% 
  mutate(id = winner) %>% 
  mutate(pre_batchcage = batch.cage)
#%>%
#  select(-batch.cage)

head(help)
tail(help)

check2 <- na.omit(help)

# 
table(help$pre_batchcage)

### add time

library(lubridate)
help<- help %>% 
  mutate(date = date(help$ts)) %>% 
  mutate(hrs = hour(help$ts)) %>% 
  mutate(min = minute(help$ts)) %>% 
  mutate(sec = second(help$ts)) %>% 
  select( -ts)  # add dates and times and remove timestamp

help <- help %>%   mutate(lightsoff = ifelse(batch %in% 1:3,11,10)) #enter lights off


help$ztime <- ((help$hrs - help$lightsoff)*3600) + (help$min*60) + help$sec  #Ztime in seconds

head(help)
tail(help)




write.csv(help, "manuscript/behavior/Pre_WL.csv", row.names=F)
