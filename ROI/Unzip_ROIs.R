###################################
# Create directories for all ROI files 
###################################
rm(list=ls())
library(oro.dicom)
library(oro.nifti)
library(plyr)
library(scales)
library(fslr)
library(cttools)
########### UNZIP ROIS.zip then run this progrma
#### delete all ROI files
### find . -regextype posix-extended -regex "^./[0-9].*[0-9]$" -exec rm -r {} \;

setup <- function(id){
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

  study = "ROI_data"


  rootdir <<- path.expand(rootdir)
  homedir <<- file.path(rootdir, study)
  homedir <<- path.expand(homedir)

#progdir <- file.path(dirname(basedir), "programs")
  progdir <<- file.path(rootdir, "programs")
  # source(file.path(progdir, "convert_DICOM.R"))
  # source(file.path(progdir, "fslhd.R"))

  basedir <<- file.path(homedir, id)

}

### initial setup
# iid <- length(ids)

#### loop through IDS and convert them to nii, gantry tilted
### 301-520 needs to use Study Date/Time instead of Series Date/Time
# for (iid in 1:length(ids)){
rerunroi = FALSE

setup("Long_ROIS")
# setup("ROIS")

    gf = getfiles(basedir)
    all.files = paths = gf$files
    bn = basename(paths)
    bnn = basename(dirname(paths))

    paths = dirname(paths)
    ids = sapply(strsplit(basename(paths), "_"), function(x) x[1])
    ids = gsub("(\\d\\d\\d)((\\d|)\\d\\d\\d)", "\\1-\\2", ids)
    iids = gsub("-", "", ids)

    new.names = paste0(bnn, "_", bn)

    id.dirs = file.path(homedir, ids)

    l_ply(unique(id.dirs), function(x){
      if (!file.exists(x)){
        dir.create(x, showWarnings=FALSE)
      }
    })

    new.names = file.path(id.dirs, new.names)
    df = data.frame(filename = all.files, 
      new.name = new.names,
      stringsAsFactors=FALSE)

    # df = df[ grepl("100-365", df$filename), ]

    file.rename(df$filename, df$new.name)
    # m_ply(function(filename, new.name){
    #   file.rename(filename, new.name)
    #   return(NULL)
    # }, .data=df, .progress = "text")

