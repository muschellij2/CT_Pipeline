###################################################################
## This code is for prediction of Image Segmentation of CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
###################################################################
###################################################################
rm(list=ls())
library(plyr)
library(cttools)
library(fslr)
library(ROCR)
library(matrixStats)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

outdir = file.path(basedir, "results")

	#### load voxel data
	outfile = file.path(outdir, "Voxel_Info.Rda")
	load(file=outfile )

	outfile = file.path(outdir, "111_Filenames_with_volumes.Rda")
	load(file = outfile)

    template.file = file.path(tempdir, "scct_unsmooth.nii.gz")
    ss.tempfile = file.path(tempdir, "Skull_Stripped",
        "scct_unsmooth_SS_First_Pass_0.1.nii.gz")
    ss.maskfile = file.path(tempdir, "Skull_Stripped",
        "scct_unsmooth_SS_First_Pass_Mask_0.1")

    tempimg = readNIfTI(ss.tempfile, reorient=FALSE)
    ssmask = readNIfTI(ss.maskfile, reorient=FALSE)


	fnames = fdf$synssimg[fdf$truevol < 20]
	imgs = t(laply(fnames, function(x){
		img = c(readNIfTI(x, reorient=FALSE))
		img = img[ssmask == 1]
		img
	}, .progress = "text"))

	blank.arr = array(NA, dim=dim(ssmask))
	mean.arr = blank.arr
	mean.arr[ssmask == 1] = rowMeans(imgs)
	mean.img = niftiarr(tempimg, mean.arr)

	sd.arr = blank.arr
	sd.arr[ssmask == 1] = rowSds(imgs)
	sd.img = niftiarr(tempimg, sd.arr)


	ortho2(mean.img)
	# x11()
	ortho2(sd.img)

	img1 = readNIfTI(fdf$synssimg[ fdf$truevol > 20][1], 
		reorient=FALSE)
	zimg1 = niftiarr(img1, (img1 - mean.img) / sd.img)
	zimg1[ is.nan(zimg1)] = NA
	zimg1[ is.infinite(zimg1)] = NA
	# x11()
	ortho2(zimg1)