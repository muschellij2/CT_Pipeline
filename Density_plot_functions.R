pos <- function(x){
  x2 <- x[x > 0]
}
scdens <- function(x, scale=TRUE){
  d <- density(x)
  if (scale) d$y <- d$y/max(d$y)
  return(d)
}
dens_plot <- function(vals, scale=TRUE, margin=NULL, showleg=TRUE, ...){
  nlines <- length(vals)
  rcols <- c("green", "yellow", "red")
  colfunc <- colorRampPalette(rcols)
  cols <- colfunc(nlines)
  dens <- lapply(vals, scdens, scale=scale)
  toplim <- max(sapply(dens, function(d) max(d$y)))
  plot(x=1:100, y=rep(0, 100), ylim = c(0, toplim), xlim=c(0, 100), type="n", xlab="Houndsfield Units", ylab=paste0("Density", ifelse(scale, ": Scaled", "")), ...)
  legs <- switch(margin,
                 "1" = c("left", "middle", "right"),
                 "2" = c("back", "middle", "front"),
                 "3" = c("bottom", "middle", "top"))
  if (showleg) legend(0, toplim, legend=legs, col=rcols, lty=1)
  for (iline in 1:nlines) lines(dens[[iline]], col=cols[iline])
}


brain_dens <- function(img, margin=3, ...){
  ### read data, take values > 0 then plot density of 
  ### margin z = 3, x= 1, and color by where in the brain
  
  
  dimg <- dim(img)
  if (length(dimg) > 3) img <- img[,,,1]
  keep <- apply(img, margin, function(x) sum(x > 0) > 1)
  
  vals <- apply(img, margin, pos)
  vals <- vals[keep]
  dens_plot(vals, margin = margin, ...)
  return(invisible(NULL))
}


