########################################################
## This code is for Image Segmentation of CT
## The code is R but calls out FSL
##
## Author: John Muschelli
## Last updated: May 20, 2014
########################################################
########################################################
rm(list=ls())
library(plyr)
library(cttools)
library(fslr)
library(mgcv)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
segdir = file.path(progdir, "Segmentation")

basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")
outdir = file.path(basedir, "results")

#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )

outfile = file.path(outdir, 
	"111_Filenames_with_volumes.Rda")
load(file = outfile)


iimg <- suppressWarnings(as.numeric(
	Sys.getenv("SGE_TASK_ID")))
if (is.na(iimg)) iimg = 34
## 2,12,17,34,37, 46, 48, 68, 70,81, 85, 86, 87, 99
#  has variable slice thickness
## 15 is not correct for sthickness
## 17 & 87 worst - has overlapping slice somewhat
## 75 has all negative
## 71 has no position data
## 13,71,101 has spacing


fdf$median = fdf$mean = fdf$mode = NA
# fdf = fdf[c(2,12,17,34,37, 46, 48, 68, 70,
	# 81, 85, 86, 87, 99),]

for (iimg in seq(nrow(fdf))){
	
	runx = x = fdf[iimg,]
	ssimg = readNIfTI(x$ssimg, reorient=FALSE)
	vals = ssimg[ssimg != 0]
	tab = table(vals)
	uvals = as.numeric(names(tab))
	fdf$median[iimg] = median(vals)
	fdf$mean[iimg] = mean(vals)
	fdf$mode[iimg] = uvals[which.max(tab)]
	# print(warnings())
	print(iimg)
}

outfile = file.path(outdir, 
	"111_Filenames_with_volumes_stats.Rda")

save(fdf, file = outfile)