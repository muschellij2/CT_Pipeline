################################
# Written 2013Oct24
# Author: John Muschelli
# Purpose: Collapse data for the 5 cases and run a multi-subject model
# Output: ROC Curve
# Use of output: prediction model for classifier
################################

rm(list=ls())
library(oro.dicom)
library(plyr)
library(AnalyzeFMRI)
library(car)
library(ROCR)
library(ggplot2)
library(data.table)
basedir <- "/Volumes/DATA_LOCAL/Image_Processing/Test_5"
if (Sys.info()[["user"]] %in% "jmuschel") basedir <- "/dexter/disk2/smart/stroke_ct/ident/Test_5"
progdir <- file.path(dirname(basedir), "programs")
source(file.path(progdir, "Zscore.R"))
source(file.path(progdir, "List_Moment_Generator.R"))
source(file.path(progdir, "pred.prob.R"))

get.volume <- function(img, thick.img){
	arr <- img * thick.img
	### don't need z-direction -- has slice thickness
	dims <- img@pixdim[2:3]
	vsize <- prod(dims)
	arr <- arr * vsize
	sum(arr)
}

files <- list.files(path=basedir, pattern="*.nii.gz", full.names=TRUE)
stubs <- gsub(".nii.gz", "", files, fixed=TRUE)
stubs <- stubs[!grepl("Zero", stubs)]
stubs <- stubs[grepl("ROI", stubs)]
mat <- matrix(stubs, ncol=1, byrow=FALSE)
mat <- data.frame(mat, stringsAsFactors=FALSE)

mat[, 2] <- gsub("_ROI", "", mat[,1])
### make sure the data only has images and ROIS
colnames(mat) <- c("roi", "img")
mat <- mat[, c("img", "roi")]
test <- gsub("_ROI", "", mat$roi)
stopifnot(all(test == mat$img))

mat$ss.img <- paste0(mat$img, "_SS_First_Pass_Mask_0.1")
mat$ss.img <- file.path(dirname(mat$ss.img), 
  "Skull_Stripped", 
  basename(mat$ss.img))

### test to see if all files exist
imgs <- unlist(mat)
imgs <- paste0(imgs, ".nii.gz")
stopifnot(all(file.exists(imgs)))


## all combinations
nfiles <- nrow(mat)
combos <- combn(nrow(mat), floor(nrow(mat)/2))
subset <- sample(1:ncol(combos), 1)
ind <- 1:nfiles

rda <- file.path(basedir, "Collaped_Data.rda")
load(file=rda)

rda <- file.path(basedir, "Collaped_Model.rda")
load(file=rda)

test.ind <- !train.ind

train <- df[samp,]
test <- df[ !(1:nrow(df) %in% samp), ]

system.time(preds <- sapply(mods, pred.prob, test=test ))
nmods <- ncol(preds)
cn <- colnames(preds) <- paste0("mod", 1:nmods)
test <- cbind(test, preds)
preds <- NULL

test <- test[, c("fname", "y", cn), with=FALSE]
test <- as.data.frame(test, stringsAsFactors=FALSE)


sums <- ddply(test, "fname", summarise,
	m1=sum(mod1),
	m2=sum(mod2), 
	# m3=sum(mod3), 
	# m4=sum(mod4), 
	# m5=sum(mod5), 
	sy=sum(y))

x <- sapply(sums[, paste0("m", 1:nmods)], function(x) (x - sums$sy))


##################################################
########## Make predicted probability maps
##########
##########
##################################################

runscans <- data.frame(img=dropscan, stringsAsFactors=FALSE)
imgs <- merge(mat, runscans, all.x=FALSE, by="img")
irow <- 1
nimgs <- nrow(imgs)
vols <- matrix(NA, nrow=nimgs, ncol=nmods)
colnames(vols) <- paste0("mod", 1:nmods)
vols <- data.frame(vols)
vols$fname <- imgs$img
vols$manual <- NA

for (irow in 1:nimgs){

	fname <- imgs$img[irow]
	folname <- basename(fname)

### get mask again (could get from test$y)
  truth <- readNIfTI(imgs$roi[irow], reorient=FALSE)
  
##read in slice thickness file
  thick.file <- file.path(basedir, "Slice_Thickness", folname)  
  thick.img <- readNIfTI(thick.file, reorient=FALSE)
  thick.img <- thick.img - thick.img@scl_inter
  vols$manual[irow] <- get.volume(truth, thick.img)

## read in brain-extracted image
  ss.mask <- readNIfTI(imgs$ss.img[irow], reorient=FALSE)
  mask <- ss.mask > 0

## 
  pimg <- ss.mask
  pimg@.Data[!mask] <- 0

  imod <- 1
  single <- test[ test$fname %in% fname, ]

  for (imod in 1:nmods){
  		modcol <- paste0("mod", imod)
 	  predictions <- single[ , modcol]
	  #all.pred <- array(all.pred, dim=dim(pimg@.Data))
	  pimg@.Data[mask] <- predictions
	  vols[irow, modcol] <- get.volume(pimg, thick.img)

	  
	  outdir <- file.path(basedir, "prob_maps", "collapsed_model")  

	  ### make a probability map (datatype=16 means float)
	  pfname <- file.path(outdir, paste0(folname, "_Model_", imod))
	  rimg <- range(pimg)
	  pimg@cal_max <- rimg[2]
	  pimg@cal_min <- rimg[1]
	  pimg@scl_inter <- 0
  	  pimg@datatype <- 16
  	  pimg@bitpix <- 32	  
	#   pimg@data_type <- "double"
	  # pimg2 <- nifti(pimg, datatype=16)
	  # pimg2@vox_offset
	  writeNIfTI(pimg, filename=pfname)
	  print(imod)
  } 
	print(irow)

}