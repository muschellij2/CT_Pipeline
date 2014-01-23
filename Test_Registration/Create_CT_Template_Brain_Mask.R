#####################################
## Author: John Muschelli
## Date: January 20, 2014
## Purpose: Read in the AAL atlas labels and make R objects
## that can be used later for overlap metrics.
#####################################
#####################################
rm(list=ls())
library(R.matlab)
library(oro.nifti)
library(plyr)
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
tempdir = file.path(rootdir, "Template")
outdir = file.path(basedir, "results")

fsldir <- system("echo $FSLDIR", intern=TRUE)
if (fsldir == "") {
	fsldir = "/usr/local/fsl"
}
spmdir = file.path(homedir, "spm8")
aaldir = file.path(spmdir, "toolbox", "aal_for_SPM8")
atlasdir = file.path(fsldir, "data", "atlases")
stddir = file.path(fsldir, "data", "standard")


template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)
tt = temp

## the r is from nii_reslice_img to get it to 1mm and correct bounding
## box
bmask.nii = file.path(stddir, "rMNI152_T1_1mm_brain_mask.nii.gz")
bmask.img = readNIfTI(bmask.nii)
stopifnot(all(bmask.img %in% c(0, 1)))
dtemp = dim(temp)
bmask.img = bmask.img > 0
temp[!bmask.img] = 0
# orthographic(temp)
t2 = temp 
t2[t2 > 100] = 0
# orthographic(t2)

tempmask = file.path(tempdir, "scct_mask")
nim = t2 > 0
nim@cal_min = min(nim, na.rm=TRUE)
nim@cal_max = max(nim, na.rm=TRUE)
nim@scl_inter = 0
nim@scl_slope = 1
writeNIfTI(nim, filename=tempmask)
## use convex hulling
fslfill(file=tempmask)
