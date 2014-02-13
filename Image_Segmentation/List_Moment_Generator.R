# rm(list=ls())
# library(oro.dicom)
# library(AnalyzeFMRI)
# homedir <- "/Users/muschellij2/Dropbox/3DPDF_Example/Paper"
# datadir <- file.path(homedir, "data")
# 
# ### Template from MNI152 from FSL
# setwd(datadir)
# template <- readNIfTI("MNI152_T1_2mm_brain.nii")

### neighbors returns list of values of neighbors, template dimensions, how the image
### was shifted for each of the list, and the indices onthe template of the data to 
### take neighbors of
neighbors <- function(image, mask=NULL, nvox=NULL){
  dimg <- dim(image) 
  if (is.null(mask)) {
    mask <- array(TRUE, dim=dimg)
  }
  if (!all.equal(dim(mask), dimg)) stop("Dimensions of mask and image do not match")
  ind <- which(mask, arr.ind=TRUE)
  image <- image * mask
  
  shift <- c(-nvox:nvox)
  shifter <- expand.grid(x=shift, y=shift, z=shift)
  ### only store x + 1, x-1, and so on once
  shifted <- lapply(shift, function(x) ind + x)
  
  ### getting the index of which to put together
  shift.ind <- t(apply(shifter, 1, match, shift))
  shift.ind <- data.frame(t(shift.ind))
  arr <- lapply(shift.ind,  function(x) {
    xcol <- x[1]
    ycol <- x[2]
    zcol <- x[3]
    cbind(x=shifted[[xcol]][,1], y=shifted[[ycol]][,2], z=shifted[[zcol]][,3])
  })
  ### make it so any edges are not included (0 or greater than FOV)
  arr <- lapply(arr, function(x){
    for (i in 1:3)  x[ x[,i] < 1 | x[,i] > dimg[i] ,i] <- NA
    x
  })
  
  vals <- lapply(arr, function(x) image[x])
  
  vals <- do.call("cbind", vals)
  return(list(vals=vals, indices=ind, dimg=dimg, shift.ind=shift.ind))
}

#### Takes the moment of an object from neighbors, and either returns the image
#### or the data scaled (somewhat useless)
neigh_moment <- function(neigh, moment=1, central=TRUE, retimg=TRUE){
  vals <- neigh$vals
  if (central & moment == 1) {
    print('Central and moment =1 will give demeaned image, central=FALSE')
    central=FALSE
  }
  
  if (central) vals <- vals - rowMeans(vals, na.rm=TRUE)
  vals <- vals^moment
  moment <- rowMeans(vals, na.rm=TRUE)
  if (!all(is.finite(moment))) stop("Someone has no neighbors?!")
  if (retimg) {
    img <- array(NA, dim=neigh$dimg)
    img[neigh$indices] <- moment
    return(img)
  } else {
    return(moment)
  }
}


### combined 2 functions to give one that allows you to do multiple moments
image_moment <- function(image, mask=NULL, nvox=NULL, moment, 
                         central=TRUE, retimg=TRUE){
  neigh <- neighbors(image, mask=mask, nvox=nvox)
  lmom <- length(moment)
  if (lmom > 1 & length(central) == 1) central <- rep(central, lmom)
  l <- list()
  for (imom in 1:lmom){
    mom <- moment[imom]
    cent <- central[imom]
    l[[imom]] <- neigh_moment(neigh, moment=mom, central=cent, retimg=retimg)
  }
  if (imom == 1) l <- l[[1]]
  return(l)
}

# 
# 
# nvox <- 1
# mask <- template > 0
# ind <- which(mask, arr.ind=TRUE)
# dind <- which(!mask, arr.ind=TRUE)
# mask[dind] <- NA
# 
# system.time(mnimg <- image_moment(image=template, mask, nvox=nvox, moment=c(1,2, 3), central=TRUE, retimg=TRUE))
# ###   user  system elapsed 
# ###   4.509   2.200   7.409
# system.time(mnimg2 <- image_moment(image=template, mask=NULL, nvox=nvox, moment=c(1,2, 3), central=TRUE, retimg=TRUE))
# ### user  system elapsed 
# ### 17.658   6.630  25.036
# 
# # 
# # neigh <- neighbors(template, mask, nvox=nvox)
# # vals <- neigh$vals
# # 
# # 
# # x <- which(neigh$indices[,1] == 54 & neigh$indices[,2] == 91 & neigh$indices[,3] == 43)
# # 
# # high <- which(template > 8330, arr.ind=TRUE)
# 




