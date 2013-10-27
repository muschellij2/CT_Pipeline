zscore <- function(mat, margin=3){
  nif_img <- FALSE
  ### check for nifti image
  if ("nifti" %in% class(mat) ) {
    xmat <- mat
    mat <- mat@.Data
    nif_img <- TRUE
  }
  
  ## scale and center the image
  for (ind in 1:dim(mat)[margin]){
    if (margin == 1) slice <- mat[ind,,]
    if (margin == 2) slice <- mat[,ind,]
    if (margin == 3) slice <- mat[,,ind]
    over0 <- slice > 0
    if (sum(over0) > 1){
      img <- slice[over0]
      mn <- mean(img, na.rm=TRUE)
      std <- sd(img, na.rm=TRUE)
      # print(mn)
      # print(std)
      slice <- (slice - mn)/std
      slice[!over0] <- 0
    }
    if (margin == 1) mat[ind,,] <- slice
    if (margin == 2) mat[,ind,] <- slice
    if (margin == 3) mat[,,ind] <- slice
  }
  if (nif_img) {
    xmat@.Data <- mat
    mat <- xmat
  }  
  mat
}


zscore2 <- function(mat, margin=3){
  require(plyr)
  ret <- aaply(.data=mat, margin, function(slice){
    over0 <- slice > 0
    if (sum(over0) > 1){
      img <- slice[over0]
      mn <- mean(img, na.rm=TRUE)
      std <- sd(img, na.rm=TRUE)
      # print(mn)
      # print(std)
      slice <- (slice - mn)/std
      slice[!over0] <- 0
    }
    slice
  })
  ret <- aperm(ret, c(2, 3, 1))
  ret
}

slice_hist <- function(mat, margin=3, plothist=TRUE, ...){
  i <- 0
  colors <- rainbow(dim(mat)[margin])
  
  for (ind in 1:dim(mat)[margin]){
    if (margin == 1) slice <- mat[ind,,]
    if (margin == 2) slice <- mat[,ind,]
    if (margin == 3) slice <- mat[,,ind]
    over0 <- slice > 0		
    if (sum(over0) > 1){
      img <- slice[over0]
      ax <- c("x", "y", "z")[margin]
      if (plothist) hist(img, main=paste0("Image intensity Density ", ax, " Slice #", ind), xlim=c(0, 100), xlab="Image Intensity of Non-Zero Voxels", breaks=seq(0, 100, 10), ...)
      else {
        d <- density(img)
        d$y <- d$y /max(d$y)
        #				color <- grey(ind/dim(mat)[margin])
        color <- colors[ind]
        if (i == 0) plot(d, main=paste0("Image intensity Density Colored by location"), xlim=c(0, 100), xlab="Image Intensity of Non-Zero Voxels", col = color)
        else lines(d, col = color)
        i = i + 1		
      }
    } else NULL	
  }
}