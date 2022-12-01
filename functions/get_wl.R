get_wl <- function(df){
  df$winner <- ifelse(df$score==0, df$Animal.2, df$Animal.1)
  df$loser <- ifelse(df$score==0, df$Animal.1, df$Animal.2)
  return(df)
}
