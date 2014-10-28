###################################################################
## This code is for Image Segmentation of CT
## The code is R but calls out FSL
##
## Author: John Muschelli
## Last updated: May 20, 2014
###################################################################
###################################################################
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

outfile = file.path(outdir, "111_Filenames.Rda")
xxx = load(file = outfile)

iimg <- suppressWarnings(as.numeric(Sys.getenv("SGE_TASK_ID")))
if (is.na(iimg)) iimg = 2
## 2 has variable slice thickness
## 15 is not correct for sthickness

fdf$posvol = fdf$varslice = fdf$gantry = fdf$truevol = NA

for (iimg in seq(nrow(fdf))){
	
	runx = x = fdf[iimg,]
	sortdir = file.path(x$iddir, "Sorted")

	# run_model = function(x, fpr.stop = .1){
	fname = xfname = nii.stub(x$img, bn=TRUE)
	rda = file.path(sortdir, paste0(fname, "_Header_Info.Rda"))
	xrda = load(rda)

	# cn = c("0018-5100-PatientPosition", 
	# 	"0020-0032-ImagePositionPatient")
	###############################
	# Slice thickness or Patient position?
	###############################
	gant = unique(as.numeric(
		dcmtables[, "0018-1120-GantryDetectorTilt"]))
	stopifnot(length(gant) == 1)

	cn = c("0020-0032-ImagePositionPatient")	
	dcmnames = colnames(dcmtables)
	stopifnot(cn %in% dcmnames)
	pos = dcmtables[, cn]
	#########################
	# 102-391 is all messed up
	#########################
	rownames(pos) = NULL
	imgpos = t(sapply(strsplit(pos, " "), as.numeric))
	colnames(imgpos) = c("x", "y", "z")
	posdiff = diff(imgpos[, "z"])
	print(posdiff)
	# posdiff = c(posdiff[1], posdiff)
	add.posdiff = c(posdiff[1], posdiff)

	tt = thicks = dcmtables[, "0018-0050-SliceThickness"]
	print(thicks)
	thicks = as.numeric(thicks)	
	if (any(is.na(thicks))){
		print(thicks)
		print(tt[is.na(thicks)])
	}
	has.varslice = length(unique(thicks)) > 1


	fname = paste0(fname, "_predictors.Rda")
	outfile = file.path(x$outdir, fname)
	xx = load(file=outfile)

	## Overlaid densities of ROIs
	# Rerun localization 
	df = img.pred$df
	stopifnot(all(df$Y %in% c(0, 1)))
	nim = img.pred$nim
	dimg = dim(nim)



	all.ind = expand.grid(lapply(dimg, seq))
	colnames(all.ind) = paste0("dim", 1:3)
	all.ind$sthick = NA

	zslices = data.frame(sthick = thicks)
	zslices$dim3 = seq(nrow(zslices))
	for (iz in seq(nrow(zslices))){
		ind = which(all.ind$dim3 == zslices$dim3[iz])
		all.ind$sthick[ind] = zslices$sthick[iz]
	}

	vdim = prod(pixdim(nim)[2:3])
	df$Y = df$Y * vdim
	sthick = all.ind$sthick
	df$sthick = sthick

	df$calcY = df$Y * df$sthick

	df$posY = df$Y * pixdim(nim)[4]

	truevol = sum(df$calcY) / 1000
	posvol = sum(df$posY) / 1000

	fname = file.path(x$outdir, 
		paste0(xfname, "_truevolume.Rda"))
	save(truevol, posvol, sthick, thicks, 
		imgpos, posdiff, add.posdiff, 
		file=fname)
	
	fdf$truevol[iimg] = truevol
	fdf$posvol[iimg] = posvol
	fdf$varslice[iimg] = has.varslice
	print(iimg)
	# print(warnings())
}

outfile = file.path(outdir, "111_Filenames_with_volumes.Rda")
save(fdf, file = outfile)
