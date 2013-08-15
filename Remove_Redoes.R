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
ids <- c("205-509", "205-517", "205-519", "225-502", "225-503", "225-504", 
"225-505", "225-506", "225-507", "225-510", "225-511", "225-515", 
"225-522", "225-523", "225-524", "232-512", "232-513", "232-514", 
"232-516", "232-521", "232-526", "289-518", "289-525", "301-520")

verbose=TRUE
untar = FALSE
convert <- TRUE
skullstrip <- TRUE
regantry <- TRUE
untgantry <- FALSE
runall <- FALSE

### initial setup
iid <- 3
id <- ids[iid]
setup(id)

for (id in ids){
  setup(id)
  sortdir <- file.path(basedir, "Sorted")
  x <- list.dirs(path=sortdir, recursive=FALSE)
  for (ifol in x){
    fname <- basename(ifol)
    if (fname %in% c("dexter", "Redo")) {
      dexfiles <- list.files(path=ifol, recursive=TRUE)
      # if (length(dexfiles) == 0) 
        system(sprintf('rm -r "%s"', ifol))
    }
  }
  print(id)
  print(x)
}
