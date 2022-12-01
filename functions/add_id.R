add_id <- function(df){
  df$id1 <- paste(df$OG.Cage.Animal.1, df$Animal.1, sep="-")
  df$id2 <- paste(df$OG.Cage.Animal.2, df$Animal.2, sep="-")
  return(df)
}