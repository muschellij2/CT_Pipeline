###################################################################
## This code is for prediciton of Image Segmentation of CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
####################################################################
####################################################################
rm(list=ls())
library(plyr)
library(cttools)
library(fslr)
library(ROCR)
library(matrixStats)
library(extrantsr)
library(ANTsR)
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

    fdf$ss = gsub("_Mask", "", fdf$mask)
    iimg <- as.numeric(Sys.getenv("SGE_TASK_ID"))
    if (is.na(iimg)) iimg = 1
    vol = fdf$truevol[iimg]


    stubfile = function(x, d = NULL, ext = "", bn=TRUE){
      b = nii.stub(x, bn=bn)
      b = paste0(b, ext)
      if (!is.null(d)){
        b = file.path(d, b)
      }
      return(b)
    }

    ####################################
    ## Run both with the Skull Stripped and not skull stripped
    ####################################

    outprefix = tempfile()

    x = fdf[iimg,]

    infile = stubfile( x$ss, x$outdir, 
        ext = "_template_zscore.nii.gz" , bn=TRUE)
    # outprefix = stubfile(x$ss, d = x$outdir)
    ofile = stubfile(infile, ext = "_smooth2", bn=FALSE)
    
    res = fslsmooth(file = infile, sigma=2, retimg = TRUE,
        mask = x$mask, outfile = ofile)    