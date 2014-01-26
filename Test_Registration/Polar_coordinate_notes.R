################################
# Written 2013Oct29
# Author: John Muschelli
# Purpose: Make NIfTI files from DICOMS
# Output: NIfTI images for images and ROIs, and slicethickness image
# Slice thickness image can be used for volume weighting
################################

rm(list=ls())
library(oro.dicom)
library(oro.nifti)
library(MASS) 
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"  
}
rootprogdir = file.path(rootdir, "programs")
## need to make package
source(file.path(rootprogdir, "fslhd.R"))
progdir = file.path(rootprogdir, "Test_Registration")
basedir = file.path(rootdir, "Test_Registration")
dirs <- list.dirs(basedir, recursive=FALSE, full.names=FALSE)
ptpat <- "\\d\\d\\d-(\\d|)\\d\\d\\d"
ids <- grep(paste0(ptpat, "$"), dirs, value=TRUE)

get.nifti.header = function(nim){
	sn = slotNames(nim)
	sn = sn[ sn != ".Data"]
	hdr = sapply(sn, slot, object=nim)
}

setwd(file.path(basedir, "RawNIfTI"))

file = "100318_20010723_0956.nii"

img = readNIfTI(file)

mask = img >= 0 & img <= 100

img[!mask] = 0
img@cal_max <- max(img, na.rm=TRUE)
img@cal_min <- min(img, na.rm=TRUE)

mask.ind = which(mask, arr.ind=TRUE)

quants = apply(mask.ind, 2, quantile, probs= c(0.2, 0.78))

brain = mask.ind
brain = brain[ 
	brain[,1] >= quants[1,1] & brain[,1] <=  quants[2,1] &
	brain[,2] >= quants[1,2] & brain[,2] <=  quants[2,2] & 
	brain[,3] >= quants[1,3] & brain[,3] <=  quants[2,3],  ]

cog1 = round(apply(brain, 2, median))
cog2 = cog1
cog2[2] = quantile(brain[,2], 0.8)

x11()
orthographic(img, xyz = cog1)
x11()
orthographic(img, xyz = cog2)

cog1
cog2
#107.800000
#-147.600000
#-45.500000



cart2pol <- function(x, y, origin=c(0, 0))
{
	x = (x- origin[1])
	y = y - origin[1]
  r <- sqrt(x^2 + y^2)
  t <- atan(y/x)

  cbind(r,t)
}

slice = mask.ind[mask.ind[,3] ==5,]
scog = colMeans(slice)[1:2]
smcog = apply(slice[, 1:2], 2, median)
slice.int = img[slice]

pol = cart2pol(slice[, 1], slice[,2], origin =scog)

crit = slice.int < 40 & slice.int > 10

select = slice[crit,]
clust = kmeans(select[, 1:2], 2)

kd = kde2d(select[,1], select[,2], 
	h = c(width.SJ(select[,1]), width.SJ(select[,1]))
	)
sscog = select[which.max(kd$z),1:2]

plot(select[ ,1], select[,2], col=clust$cluster)
points(scog[1], scog[2], col="red", pch=16, cex=2)
# points(smcog[1], smcog[2], col="green", pch=16, cex=2)
points(sscog[1], sscog[2], col="green", pch=16, cex=2)


plot(pol[, "r"], y= slice.int, pch='.')
keep = pol[, "r"] < 145

plot(slice[keep, 1:2])
points(scog[1], scog[2], col="red", pch=16, cex=2)