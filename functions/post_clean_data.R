post_clean_data <- function(df){
  df <- df[!is.na(df$score),]
  df <- df[df$Start.End!="Start",]
  df <- df[df$Start.End!="End",]
  df <- df[df$id1!=df$id2,]
  
}