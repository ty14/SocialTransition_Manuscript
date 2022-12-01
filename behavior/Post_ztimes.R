library(tidyverse)


# source
df <- read.csv("data_raw/Post_startTimes.csv")
head(df)


dfl <- df %>% pivot_longer(cols = 2:3, names_to = "times")
head(dfl)

## extracting data 
library(lubridate)
dfl$value<- as.POSIXlt(dfl$value, tryFormats = c("%H:%M:%OS"))

dfz <-dfl %>% 
  mutate(hrs = hour(dfl$value)) %>% 
  mutate(min = minute(dfl$value)) %>% 
  mutate(sec = second(dfl$value)) %>% 
  select( -value)

head(dfz)


##Making ztimes 
dfz<- dfz %>%   mutate(lightsoff = ifelse(batch %in% 1:3,11,10)) 

dfz$ztime <- ((dfz$hrs - dfz$lightsoff)*3600) + (dfz$min*60) + dfz$sec  #Ztime in seconds
head(dfz)

write.csv(dfz, "manuscript/behavior/post_startTimes.csv",row.names=F)



