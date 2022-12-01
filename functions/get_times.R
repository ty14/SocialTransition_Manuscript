get_times <- function(dataf){
  dataf$ts <- strptime(as.character(dataf$Timestamp),'%m/%d/%Y %H:%M:%S')
  dataf <- dataf[order(dataf$ts),]
  return(dataf)
}