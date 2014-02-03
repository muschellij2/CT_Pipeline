rm(list=ls())
library(oro.dicom)
library(oro.nifti)
library(plyr)
library(scales)
library(reshape2)

#### delete all ROI files
### find . -regextype posix-extended -regex "^./[0-9].*[0-9]$"
###  -exec rm -r {} \;

  username <- Sys.info()["user"][[1]]

  cluster=FALSE
  if (username == "muschellij2"){
    # rootdir <- "/Volumes/DATA/New_Age_Test"
    rootdir <- "~/CT_Registration"
  } else {
    rootdir <- "/dexter/disk2/smart/stroke_ct/ident"
    cluster =TRUE;
  }
    rootdir <- path.expand(rootdir)

  # ss <- as.numeric(strsplit(id, "-")[[1]][2])
  # if (ss > 4000){
  #   study <- "CLEAR_III"
  #   dpath <- file.path("CLEAR", "CLEAR III")
  # } else if (ss > 300 & ss < 500){
  #   dpath <- study <- "MISTIE"
  # } else if (ss > 500 & ss < 4000) {
  #   dpath <- study <- "ICES" 
  # }


  rootdir <<- path.expand(rootdir)

#progdir <- file.path(dirname(basedir), "programs")
  progdir <- file.path(rootdir, "programs")
  source(file.path(progdir, "convert_DICOM.R"))
  source(file.path(progdir, "fslhd.R"))


#### setting up if things are on the cluster or not

load(file.path(rootdir, "Registration", 
  "Registration_Image_Names.Rda"))

imgs = mlply(.fun = function(outfile, roi.nii, raw, ss){
  c(outfile, roi.nii, raw, ss)
}, .data=df[, c("outfile", "roi.nii", "raw", "ss")])
attr(imgs, "split_labels") <- NULL

x = imgs[[1]]
rets = laply(.data=imgs, .fun = function(x){
  acpc_reorient(infiles = x)
})