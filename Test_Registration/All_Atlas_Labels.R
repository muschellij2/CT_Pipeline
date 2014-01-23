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
progdir = file.path(rootdir, "programs", "Test_Registration")
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

bmask.nii = file.path(stddir, "rMNI152_T1_1mm_brain_mask.nii.gz")
bmask.img = readNIfTI(bmask.nii)
stopifnot(all(bmask.img %in% c(0, 1)))
dtemp = dim(temp)
bmask.img = bmask.img > 0
temp[!bmask.img] = 0
orthographic(temp)
t2 = temp 
t2[t2 > 100] = 0
orthographic(t2)

tempmask = file.path(tempdir, "scct_mask")
nim = t2 > 0
nim@cal_min = min(nim, na.rm=TRUE)
nim@cal_max = max(nim, na.rm=TRUE)
nim@scl_inter = 0
nim@scl_slope = 1
writeNIfTI(nim, filename=tempmask)

# bm = bmask.img
# bmask.img = bmask.img[1:dtemp[1], 1:dtemp[2], 1:dtemp[3]]
# bmask.img = bmask.img > 0
# temp[!bmask.img] = 0

# t2 = tt
# bm = bm[2:(dtemp[1]+1), 2:(dtemp[2]+1), 2:(dtemp[3]+1)]
# bm = bm > 0
# t2[!bm] = 0

# orthographic(t2)
