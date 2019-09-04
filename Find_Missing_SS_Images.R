rm(list=ls())
library(cttools)
library(fslr)
library(plyr)
library(methods)
options(matlab.path='/Applications/MATLAB_R2014b.app/bin')

setup <- function(id, study = "Registration"){
  username <- Sys.info()["user"][[1]]

  cluster=FALSE
  if (username %in% c("muschellij2", "johnmuschelli")){
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
  homedir <<- file.path(rootdir, study)
  homedir <<- path.expand(homedir)

#progdir <- file.path(dirname(basedir), "programs")
  progdir <- file.path(rootdir, "programs")
  # source(file.path(progdir, "convert_DICOM.R"))
  # source(file.path(progdir, "fslhd.R"))

  basedir <<- file.path(homedir, id)

}


#### setting up if things are on the cluster or not
## ROIformat after 134-327.zip
ROIformat = TRUE
study = "Registration"
if (ROIformat) {
  study = "ROI_data"
}

setup(study, study=study)

if (ROIformat) setup("Long_ROIS", study=study)


ids = list.dirs(homedir, recursive=FALSE, full.names=FALSE)
ids = basename(ids)
ids = grep("\\d\\d\\d-(\\d|)\\d\\d\\d", ids, value=TRUE)
length(ids)



iid <- as.numeric(Sys.getenv("SGE_TASK_ID"))

# iid = 100 needs ss

if (is.na(iid)) iid <- 118


all.df = NULL
# for (iid in seq_along(ids)){
for (iid in seq_along(ids)){
  
  id <- ids[iid]
  print(id)
  setup(id, study = study)
  # source(file.path(progdir, "file_functions.R"))


  iddir = file.path(rootdir, "Registration", id)
  ssdir = file.path(iddir, "Skull_Stripped")

  imgs = list.files(iddir, pattern=".nii.gz", recursive=FALSE)

  ssimgs = file.path(ssdir, gsub("[.]nii", "_SS_0.01_Mask.nii", imgs))

  df = data.frame

  df = data.frame(img=imgs, ssimg = ssimgs, stringsAsFactors=FALSE)
  all.df = rbind(all.df, df)  
}

all.df = all.df[ !file.exists(all.df$ssimg),]
