

Heri.asreml <- function(Model,genotype){
  
  # Genetic variance component
  vc.g <- summary(Model)$varcomp[paste0(genotype,"!",genotype,".var"),'component']
  vc.g #0.142902
  
  # Mean variance of a difference of two genotypic BLUPs
  vdBLUP.mat <- predict(Model, classify=genotype, only=genotype, sed=TRUE)$pred$sed^2 # obtain squared s.e.d. matrix 
  vdBLUP.avg <- mean(vdBLUP.mat[upper.tri(vdBLUP.mat, diag=FALSE)]) # take mean of upper triangle
  vdBLUP.avg 
  
  #############
  # H2 Cullis #
  #############
  H2Cullis <- 1 - (vdBLUP.avg / 2 / vc.g)
  
  return(H2Cullis)
  
}


# For spatial design 
Heri.AR.AR <- function(Model,Genotype){
  
  sp4.pv.h2 <- predict(Model, classify = Genotype, only = Genotype, sed = TRUE, control = list(trace = FALSE))
  avepev <- mean(sp4.pv.h2$predictions$pvals$standard.error^2)
  vcV <- sp4.pv.h2$gammas["line!line.var"]*sp4.pv.h2$sigma2
  h2.wn4 <- 1 - avepev/(vcV); h2.wn4[1]
  
  return(h2.wn4)
}
