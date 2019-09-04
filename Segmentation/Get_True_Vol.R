
#####################################################
## This code is for Image Segmentation of CT
## The code is R but calls out FSL
##
## Author: John Muschelli
## Last updated: May 20, 2014

#####################################################
#####################################################
rm(list=ls())
library(plyr)
library(cttools)
library(fslr)
library(mgcv)

library(methods)
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


iimg <- suppressWarnings({
	as.numeric(Sys.getenv("SGE_TASK_ID"))
	})
if (is.na(iimg)) iimg = 34
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

# fdf = fdf[c(2,12,17,34,37, 46, 48, 68, 70,81, 
	#85, 86, 87, 99),]
# dcmtables[, '0018-1152-Exposure']


for (iimg in seq(nrow(fdf))){
	
	runx = x = fdf[iimg,]
	sortdir = file.path(x$iddir, "Sorted")

	# run_model = function(x, fpr.stop = .1){
	fname = xfname = nii.stub(x$img, bn=TRUE)

	rda = file.path(sortdir, paste0(fname, 
		"_Header_Info.Rda"))

	xrda = load(rda)
	# print(grep("pac", colnames(dcmtables), value=TRUE))
	# print(iimg)

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
	print(unique(grep("xposure", dcmnames, value=TRUE)))


	stopifnot(cn %in% dcmnames)
	pos = dcmtables[, cn]
	if (any(pos == "")) {
		cat(paste0(x$id, " has no position data\n"))
	}
	pos[ pos == ""] = "NA NA NA"



	#########################
	# 102-391 is all messed up
	#########################
	imgpos = t(sapply(strsplit(pos, " "), as.numeric))
	colnames(imgpos) = c("x", "y", "z")
	dimgpos = diff(imgpos)

	pz = imgpos[, "z"]
	# ord = order(pz)
	# dcmtables = dcmtables[ord,]
	# pz = pz[ord]
	# imgpos = imgpos[ord, ]
	####################

	# Must use pythag theorem for oblique acquisition 
	# not just diff
	### but should use diff if orthogonal - for example #
	# if negative
	# z direction means overlapping slice - 
	# shouldnt' count twice

	####################
	posdiff = sqrt(rowSums(dimgpos^2))
	tester = c(dimgpos[,1:2])
	if (all(abs(tester) <= 1e-6 & !is.na(tester))){
		posdiff = diff(pz)
	}

	if (all(posdiff < 0 & !is.na(posdiff))) {
		posdiff = -posdiff
	}
	# print(posdiff )
	add.posdiff = c(posdiff, posdiff[length(posdiff)])


	orient = dcmtables[, 
	"0020-0037-ImageOrientationPatient"]



	tt = thicks = dcmtables[, "0018-0050-SliceThickness"]
	# print(thicks)
	thicks = as.numeric(thicks)	
	if (any(is.na(thicks))){
		print(thicks)
		print(tt[is.na(thicks)])
	}
	has.varslice = length(unique(thicks)) > 1

	ind = which(diff(thicks) != 0)
	thick.diff= which(abs(thicks - add.posdiff) > 1)
	ind = sort(unique(c(ind, thick.diff)))
	if (length(ind) > 0 ){

		ind = c(max(1, ind-1), ind, min(ind+1, 
			length(thicks)))

		print(cbind(thicks, add.posdiff)[ind,])
		print(orient[1])
		print(iimg)
	}
	if (all(is.na(add.posdiff))){
		add.posdiff = thicks
	}
# }	
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
	all.ind$pos = all.ind$sthick = NA

	zslices = data.frame(sthick = thicks, pos=add.posdiff)
	zslices$dim3 = seq(nrow(zslices))
	for (iz in seq(nrow(zslices))){
		ind = which(all.ind$dim3 == zslices$dim3[iz])
		all.ind$sthick[ind] = zslices$sthick[iz]
		all.ind$pos[ind] = zslices$pos[iz]
	}

	###############################
	# Calculate like osirix
	###############################
	df$dim3 = all.ind$dim3
	vols = ddply(df, .(dim3), summarise, vol=sum(Y))
	# rep.val = 0
	rep.val = NA
	vols$vol[vols$vol == 0] = rep.val
	vols$v = c(rep.val, vols$vol[seq(nrow(vols)-1)])
	vols$avg = (vols$vol + vols$v)/2
	vols = merge(vols, zslices, by="dim3")

	vdim = prod(pixdim(nim)[2:3])

	vols$avg = vols$avg * vdim 
	vols$tavg = vols$avg * vols$sthick
	vols$pavg = vols$avg * vols$pos

	ovols = colSums(vols[, c("tavg", "pavg")], 
		na.rm=TRUE)/1000


	df$Y = df$Y * vdim
	
	sthick = all.ind$sthick
	df$sthick = sthick

	allpos = all.ind$pos
	df$pos = allpos

	df$thickY = df$Y * df$sthick

	df$posY = df$Y * df$pos

	df$zY = df$Y * pixdim(nim)[4]

	truevol = sum(df$posY) / 1000
	thickvol = sum(df$thickY) / 1000
	zvol = sum(df$zY) / 1000

	fname = file.path(x$outdir, 
		paste0(xfname, "_truevolume.Rda"))
	save(truevol, thickvol, zvol, sthick, allpos, 
		thicks, 
		imgpos, posdiff, add.posdiff, 
		orient,
		ovols,
		file=fname)
	
	fdf$truevol[iimg] = truevol
	fdf$thickvol[iimg] = thickvol	
	fdf$zvol[iimg] = zvol
	fdf$varslice[iimg] = has.varslice
	fdf$gantry[iimg] = gant
	print(iimg)
	# print(warnings())
}


outfile = file.path(outdir, 
	"111_Filenames_with_volumes.Rda")

save(fdf, file = outfile)
