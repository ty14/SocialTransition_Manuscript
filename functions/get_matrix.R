get_matrix <- function(df){
  df <-df[df$score==1,]
  mat <-compete::org_matrix(compete::get_wl_matrix(df[,c('winner','loser')]),
                            method = "ds")
  return(mat)
}