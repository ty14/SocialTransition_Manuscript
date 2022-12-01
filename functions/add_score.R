add_score <- function(df){
  df$behav <-  paste(df$Animal.1..Behavior, df$Animal.2..Behavior, sep=" --- ")
  df$score <- behavs$X[match(df$behav, behavs$behavior)]
  return(df)
}
