###################################################################
## This code is for Image Segmentation of CT
## The code is R but calls out FSL
##
## Author: John Muschelli
## Last updated: May 20, 2014
###################################################################
###################################################################
rm(list=ls())
library(fslr)
library(drammsr)

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

outfile = file.path(outdir, "111_Filenames.Rda")
xxx = load(file = outfile)

iimg <- suppressWarnings(as.numeric(Sys.getenv("SGE_TASK_ID")))
if (is.na(iimg)) iimg = 6
## 2,12,17,34,37, 46, 48, 68, 70,81, 85, 86, 87, 99
#  has variable slice thickness
## 15 is not correct for sthickness
## 17 & 87 worst - has overlapping slice somewhat
## 75 has all negative
## 71 has no position data
## 13,71,101 has spacing

fdf$median = fdf$mode = fdf$mean = NA

fdf$thickvol = fdf$zvol = fdf$varslice = 
fdf$gantry = fdf$truevol = NA

# fdf = fdf[c(2,12,17,34,37, 46, 48, 68, 70,81, 85, 86, 87, 99),]
# dcmtables[, '0018-1152-Exposure']
exposures = vector(mode="list", length=nrow(fdf))

outdef = tempfile(fileext ='.nii.gz')

# for (iimg in seq(nrow(fdf))){
	tdir = "/dexter/disk2/smart/stroke_ct/ident/Template"
	fname = fdf$ssimg[iimg]
	outname = paste0(nii.stub(fname), "_ventricle_warp")

	d = dramms(target=fname,
	source=file.path(tdir, "Skull_Stripped", 
		"scct_unsmooth_SS_0.01.nii.gz"), 
	outdef = outdef, 
	retimg = TRUE)

	vent = dramms_warp(def = outdef, 
	source=file.path(tdir, "scct_unsmooth_ventricle_mask.nii.gz"), 
	retimg=TRUE)

	vent = drop_img_dim(vent)
	print(dim(vent))
	d = dim_(vent)
	dim_(vent)[d < 1] = 1

	pdim = pixdim(vent)
	pixdim(vent)[pdim < 1] = 1

	print(dim_(vent))
	print(pixdim(vent))
	writeNIfTI(vent, filename = outname)	
# }

