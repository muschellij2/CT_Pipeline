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
library(fslr)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
outdir = file.path(basedir, "results")
get.fsl()
fsldir = Sys.getenv("FSLDIR")
fsltemp = file.path(fsldir, "data", "standard")
whichdir = "reoriented"


spmdir = file.path(homedir, "spm8")

atlasdir = file.path(tempdir, "atlases")

outdim = c(181, 217, 181)




### contents
# tal.df, tal.img, 
#   mni.df, mni.img, 
#   hoxcort.df, hoxcort.img,
#   hoxsubcort.df, hoxsubcort.img, 

