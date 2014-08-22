rm(list=ls())
library(cttools)
library(fslr)
library(plyr)
homedir = "~/"
rootdir = "~/CT_Registration/"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
homedir = file.path(rootdir, "Rorden_data")
progdir = file.path(rootdir, "programs")
basedir = file.path(rootdir, "Registration")

# scandir = file.path(homedir, "ct_scans")
regdir = file.path(homedir, "registered_ct_scans")
lesdir = file.path(homedir, "CT_lesions")
# setwd(homedir)
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")
outdir = file.path(basedir, "results")

template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)
dtemp = dim(temp)


maskfile = file.path(tempdir, "scct_mask.nii.gz")
tempmask = readNIfTI(maskfile)


fname = file.path(regdir, "Mean_Image.nii.gz")
mn.img = readNIfTI( fname)

fname = file.path(regdir, "SD_Image")
sd.img = readNIfTI(fname)

fname = file.path(regdir, "CH_SD_Image")
ch_sd.img = readNIfTI(fname)

fname = file.path(regdir, "N_Image")
N.img = readNIfTI(fname)


files = list.files(path=lesdir, full.names=TRUE, 
	pattern="^w.*_CT\\.nii")

ss_files = gsub("_CT\\.", "_CT_SS_0.01_Mask.", files)

voi_files = gsub("_CT\\.", "_voi.", files)


maskimg = function(img, mask, replacer = NA){
	img[mask == 0] = replacer
	img
}


zscore = function(img, ssimg, cutoff = 10){

	img = check_nifti(img)
	ssimg = check_nifti(ssimg)
	ssimg = cal_img(ssimg > .5)

	dimg = (img - mn.img)
	zarr = dimg /sd.img
	ch_zarr = dimg/ch_sd.img

	zarr[ which(zarr > cutoff)] = cutoff
	zarr[ which(zarr < -cutoff)] = -cutoff

	ch_zarr[ which(ch_zarr > cutoff)] = cutoff
	ch_zarr[ which(ch_zarr < -cutoff)] = -cutoff

	zimg = temp
	zimg@.Data = zarr
	zimg = cal_img(zimg)

	ch_zimg = temp
	ch_zimg@.Data = ch_zarr
	ch_zimg = cal_img(ch_zimg)

	ch_zimg = maskimg(ch_zimg, ssimg)
	zimg = maskimg(zimg, ssimg)
	return(list(zimg=zimg, ch_zimg = ch_zimg))
}


ifile = 1
ssimg = ss_files[ifile]
img = files[ifile]
voiimg = voi_files[ifile]
cutoff = 10

voiimg = check_nifti(voiimg)

zs = zscore(img, ssimg)

