rm(list=ls())
library(oro.dicom)
library(bitops)
library(arules)

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

  ss <- as.numeric(strsplit(id, "-")[[1]][2])
  if (ss > 4000){
    study <- "CLEAR_III"
    dpath <- file.path("CLEAR", "CLEAR III")
  } else if (ss > 300 & ss < 500){
    dpath <- study <- "MISTIE"
  } else if (ss > 500 & ss < 4000) {
    dpath <- study <- "ICES" 
  }


  rootdir <<- path.expand(rootdir)
  homedir <<- file.path(rootdir, study)
  homedir <<- path.expand(homedir)

#progdir <- file.path(dirname(basedir), "programs")
  progdir <<- file.path(rootdir, "programs")
  source(file.path(progdir, "convert_DICOM.R"))

  basedir <<- file.path(homedir, id)

}


#### setting up if things are on the cluster or not
# id <- "238-4136"
id <- "301-520"
setup(id)
# source(file.path(progdir, "file_functions.R"))


files <- dir(path=basedir, pattern="_ungantry.tar.gz$", full.names=TRUE, 
      recursive=TRUE)
### untarball the files using R's command
if (length(files) > 0){
  for (ifile in files){
    dname <- dirname(ifile)
    untar(ifile, exdir = dname, compressed='gzip')
    if (verbose) print(ifile)
  }
}
