## Getting weight and everything together for gene analysis need continuous and cateorgical data 
library(tidyverse)

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

data <- na.omit(data)

write.csv(data, "manuscript/cort/physiologyData.csv", row.names = F)
