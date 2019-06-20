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
