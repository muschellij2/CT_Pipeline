
rm(list=ls())
library(oro.dicom)
library(oro.nifti)
library(plyr)

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

  study = "Registration_ICES"
  todir <- file.path(rootdir, study)
  fromstudy = "ICES_raw"
  fromdir = file.path(rootdir, fromstudy)

  files = list.files(fromdir, pattern="*.zip", full.names=TRUE)
  ifile = 1
  for (ifile in seq_along(files)){
    print(ifile)
    unzip(files[ifile], exdir = todir, overwrite=TRUE)
  }