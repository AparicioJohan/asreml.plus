

# Pseudo R-square
r.asreml <- function(Model,data,response){
  
  response      <- data[,response]
  mean.response <- mean(response,na.rm = T)
  fitted        <- fitted(Model)
  SS.fitted     <- sum( (response-fitted)^2 , na.rm = T)
  SS.response   <- sum( (response-mean.response)^2 , na.rm = T )
  
  R <- 1- SS.fitted/SS.response
  
  names(R) <- "r.square"
  return(round(R,3))
}