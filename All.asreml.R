library(asreml)


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



# Heritability classical design
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


# Heritability for spatial design 
Heri.AR.AR <- function(Model,Genotype){
  
  sp4.pv.h2 <- predict(Model, classify = Genotype, only = Genotype, sed = TRUE, control = list(trace = FALSE))
  avepev <- mean(sp4.pv.h2$predictions$pvals$standard.error^2)
  vcV <- sp4.pv.h2$gammas["line!line.var"]*sp4.pv.h2$sigma2
  h2.wn4 <- 1 - avepev/(vcV); h2.wn4[1]
  
  return(h2.wn4)
}



# spatial plot with nugget effect

plot.spatial.asreml <- 
  function(x,Datos,res="YdHa", col="col",row="row",gen="line",all.in.one = TRUE, main = NULL, annotated = FALSE, depict.missing = FALSE, ...) {
  
  xlab <- col
  ylab <- row
  x.coord <- Datos[,col]
  y.coord <- Datos[,row]
  response <- Datos[,res] 
  Datos$rowname <- rownames(Datos)
  
  residuals <- residuals(x,spatial = "plot")
  
  fitted <-  x$fitted.values
  
  
  geno.pred <- coef(x)$random %>% data.frame(gen=rownames(.),BLUP=.)
  geno.pred <- geno.pred[grep(pattern = gen,geno.pred$gen),]
  geno.pred$gen <- str_remove(rownames(geno.pred), paste0(gen,"_"))
  
  names(geno.pred)[1] <- gen
  
  columns <- seq(min(x.coord), max(x.coord), by = min(diff(sort(unique(x.coord)))))
  rows <- seq(min(y.coord), max(y.coord), by = min(diff(sort(unique(y.coord)))))
  xy.coord <- data.table(expand.grid(columns = columns, rows = rows))
  
  v <- merge(geno.pred,Datos,by = gen,all.y = TRUE)
  xx <- arrange(v,col,row)
  geno.pred <- as.vector(xx$effect)
  
  environment <- fitted-coef(x)$fixed["(Intercept)",] - geno.pred - coef(x, pattern = "units")[,1]   
  
  setNumericRounding(2)
  
  if(is.null(main)) main = paste("Trait: ",res , sep = "")
  
  setkeyv(xy.coord, c("rows", "columns"))
  ONE <- rep(1, length(x.coord))    
  df <- data.table(columns = x.coord, rows = y.coord, 
                   response = response, fitted = fitted,environment,
                   residuals = residuals, geno.pred = geno.pred, ONE = ONE)
  setkeyv(df, c("rows", "columns"))
  df <- df[xy.coord]
  df <- df[order(df$columns, df$rows),]
  
  colors = topo.colors(100)
  
  main.legends <- c('Raw data', 'Fitted data', 'Residuals',"Effect Design"  ,"Genotypic BLUPs", 'Histogram')
  if(all.in.one) {
    op <- par(mfrow = c(2,3), oma = c(ifelse(annotated, 12, 2), 1, 3, 2), mar = c(2.7, 4, 2.5, 2.5), mgp = c(1.7, 0.5, 0))                
  } else {
    if(!is.null(main))
      main.legends <- rep(main, length(main.legends))
  }
  
  range <- range(c(response, fitted), na.rm = TRUE)
  fields::image.plot(columns, rows, t(matrix(df$response, ncol = length(columns), nrow = length(rows))), main = main.legends[1], col = colors, xlab = xlab, ylab = ylab, zlim = range, graphics.reset = TRUE, ...)
  if(!all.in.one)
    readline("Press return for next page....")
  fields::image.plot(columns, rows, t(matrix(df$fitted, ncol = length(columns), nrow = length(rows))), main = main.legends[2], col = colors, xlab = xlab, ylab = ylab, zlim = range, graphics.reset = TRUE, ...)
  if(!all.in.one)
    readline("Press return for next page....")
  fields::image.plot(columns, rows, t(matrix(df$residuals, ncol = length(columns), nrow = length(rows))), main = main.legends[3], col = colors, xlab = xlab, ylab = ylab, graphics.reset = TRUE, ...)
  if(!all.in.one)
    readline("Press return for next page....")
  fields::image.plot(columns, rows, t(matrix(df$environment, ncol = length(columns), nrow = length(rows))), main = main.legends[4], col = colors, xlab = xlab, ylab = ylab, graphics.reset = TRUE, ...)
  if(!all.in.one)
    readline("Press return for next page....")
  fields::image.plot(columns, rows, t(matrix(df$geno.pred, ncol = length(columns), nrow = length(rows))), main = main.legends[5], col = colors, xlab = xlab, ylab = ylab, graphics.reset = TRUE, ...)
  if(!all.in.one)
    readline("Press return for next page....")
  # 
  suppressWarnings(hist(unique(geno.pred), main = main.legends[6], xlab = main.legends[6], ...))        
  title("")
  mtext(main, cex = 1.5, outer = TRUE, side = 3)
  invisible(df)
  
}
