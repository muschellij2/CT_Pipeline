rm(list=ls())
library(cttools)
library(fslr)
library(plyr)
options(matlab.path='/Applications/MATLAB_R2016a.app/bin')


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

study = "Registration"

setup(study, study=study)


ids = list.dirs(homedir, recursive=FALSE, full.names=FALSE)
ids = basename(ids)
ids = grep("\\d\\d\\d-(\\d|)\\d\\d\\d", ids, value=TRUE)
length(ids)


iid <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if (is.na(iid)) iid <- 121

id <- ids[iid]
setup(id, study = study)


# for (iid in seq_along(ids)){

  id <- ids[iid]
  print(id)
  setup(id, study = study)
  # source(file.path(progdir, "file_functions.R"))
  imgs = list.files(path=basedir, pattern = ".nii.gz", 
  	recursive = FALSE,
  	full.names=TRUE)

  imgs = imgs[!is.na(imgs)]
  print(length(imgs))
  for (iimg in seq_along(imgs)){
  	img = readNIfTI(imgs[iimg], reorient=FALSE)
  	cat(paste0(imgs[iimg], "\n"))
  	print(iimg)
  	ortho2(img)
  	# image(img, plot.type="single", plane = "sagittal", z = 200)
  	imgs[iimg] = NA
  	# Sys.sleep(2)
  }

  dev.off()