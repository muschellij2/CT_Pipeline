#################################
# Regressions with % of ROI
# Author: John Muschelli
#################################
rm(list=ls())
library(oro.nifti)
library(plyr)
library(scales)
library(RColorBrewer)
library(data.table)
library(cttools)
library(fslr)
library(ggplot2)
library(grid)
library(data.table)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
# source(file.path(progdir, "convert_DICOM.R"))
# source(file.path(progdir, "fslhd.R"))
basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

outdir = file.path(basedir, "results")


load(file=file.path(rootdir, "Registration", 
  "Registration_Image_Names.Rda"))

df = df[, c("raw", "roi.nii", "ss", "copydir")]
df$raw = paste0(df$raw, ".gz")
df$roi.nii = paste0(df$roi.nii, ".gz")
df$ss = paste0(df$ss, ".gz")
irow = 4

df$dists = NA
# winds <- -1:1	
# indices <- permutations(length(winds), 3, v = winds, 
# 	repeats.allowed=TRUE)
# indices = indices[!apply(indices, 1, function(x) all(x == 0)),]

for (irow in seq(nrow(df))){

	# raw = readNIfTI(df$raw[irow], reorient=FALSE)
	img = readNIfTI(df$raw[irow], reorient=FALSE)
	ss = readNIfTI(df$ss[irow], reorient=FALSE)
	ero = fslerode(df$ss[irow], kopts = "-kernel box 1", 
		retimg=TRUE)
	diff = ss - ero
	diff.img = ss
	diff.img@.Data = diff
	diff.img = cal_img(diff.img)

	stopifnot( all(diff.img %in% c(0, 1)))
	ind = which(diff.img > 0, arr.ind=TRUE)

	roi = readNIfTI(df$roi.nii[irow], reorient=FALSE)
	ero.roi = fslerode(df$roi.nii[irow], kopts = "-kernel box 1", 
		retimg=TRUE)
	diff.roi = diff.img
	diff.roi@.Data = roi - ero.roi

	allroi.ind = roi.ind = which(diff.roi > 0, arr.ind=TRUE)

	vdim = voxdim(roi)
	#### row-wise matrix mulitplication
	ind = t(t(ind)*vdim)
	roi.ind = t(t(roi.ind)*vdim)

	ind2 = rowSums(ind^2)
	roi.ind2 = rowSums(roi.ind^2)
	mid = -2* tcrossprod(ind, roi.ind)
	res = t(t(mid) + roi.ind2)
	res = res + ind2

	# checki = sample(nrow(ind), size=1)
	# checkj = sample(nrow(roi.ind), size=1)
	# check = res[checki, checkj]
	# r = roi.ind[checkj,]
	# i = ind[checki,]
	# d = sum((r-i)^2)
	# stopifnot(abs(d - check) < 1e-10)

	res = sqrt(res)

	### res is in mm, convert to cm
	res = res / 10
	mres = min(res)
	min.ind = which(abs(res- mres) < 1e-10, arr.ind=TRUE)

	roi[min.ind] = 100
	index = allroi.ind[min.ind[,2], ]
	index = index[1,]

	stub = nii.stub(basename(df$raw[irow]))
	pngname = file.path(df$copydir[irow], "results", 
		paste0(stub, "_min_dist.png"))
	print(paste0("Pngname is "))
	print(pngname)
	# png(pngname)
	# mask.overlay(img, roi, col.y=c("blue", "red"), 
	# 	xyz= round(index))
	# dev.off()

	# ind = data.frame(which(ss > 0, arr.ind=TRUE))
	# l = lapply(indices, function(x) x + ind)
	# ind.df = rbindlist(l)
	# setkey(ind.df)
	# print(nrow(ind.df))
	# unique(ind.df)
	# ind.df = unique(ind.df)

	# dind = data.table(ind)
	# setkey(dind)
	# dind[, inbrain := 1]
	# d = merge(ind.df, dind, all=TRUE)

	# d = d[ which(is.na(inbrain)), ]

	# mat = ss
	# mat@.Data[is.na(mat) | !is.na(mat)] = NA
	# dd = as.matrix(d)[, c("dim1", "dim2", "dim3")]
	# mat@.Data[dd] = 1
	# mat = cal_img(mat)

	df$dists[irow] = mres
	print(irow)
	# mask.overlay(ss, mat)
}

save(df, file=file.path(rootdir, "Registration", 
  "Minimum_Distance_To_Cortex.Rda"))