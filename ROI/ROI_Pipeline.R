rm(list=ls())
library(oro.dicom)
library(oro.nifti)
library(plyr)
library(scales)

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
  source(file.path(progdir, "convert_DICOM.R"))
  source(file.path(progdir, "fslhd.R"))

  basedir <<- file.path(homedir, id)

}

#### setting up if things are on the cluster or not
verbose =TRUE
untar = TRUE
convert <- TRUE
skullstrip <- FALSE
plotss = TRUE
regantry <- FALSE
untgantry <- FALSE
runall <- TRUE
useRdcmsort= TRUE
useRdcm2nii= FALSE
removeDups = TRUE
ROIformat = TRUE
dcm2niicmd = "dcm2nii_2009"

### initial setup
# iid <- length(ids)

#### loop through IDS and convert them to nii, gantry tilted
### 301-520 needs to use Study Date/Time instead of Series Date/Time
# for (iid in 1:length(ids)){
setup("ROIS")

ids = list.dirs(homedir, recursive=FALSE, full.names=FALSE)
ids = grep("\\d\\d\\d-(\\d|)\\d\\d\\d", ids, value=TRUE)
length(ids)

iid <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if (is.na(iid)) iid <- 4


for (iid in seq_along(ids)){

  id = ids[iid]

  setup(id)


  if (convert) {
    ### convert the dicoms
    infofile <- file.path(basedir, "Dropout_Information.Rda")
    file.remove(infofile)
    
    ### time for convertsion
    contime <- NULL



    contime <- system.time(convert_DICOM(basedir, progdir, 
                            verbose=verbose, untar=untar, 
                            useRdcmsort= useRdcmsort, 
                            useRdcm2nii= useRdcm2nii,
                            id = id, 
                            removeDups=removeDups,
                            ROIformat = ROIformat, 
                            dcmsortopt=dcmsortopt, 
                            dcm2niicmd=dcm2niicmd))

    print(contime)

  }


}## end for loop

