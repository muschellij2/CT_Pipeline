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


#### setting up if things are on the cluster or not
ROIformat = FALSE
study = "MISTIE_MRI"
if (ROIformat) {
  study = "ROI_data"
}

setup(study, study=study)

if (ROIformat) setup("Long_ROIS", study=study)


ids = list.dirs(homedir, recursive=FALSE, full.names=FALSE)
ids = basename(ids)
ids = grep("\\d\\d\\d-(\\d|)\\d\\d\\d", ids, value=TRUE)
length(ids)


runonlybad = FALSE

if (runonlybad) ids = ids[bad.ids]

verbose =TRUE
untar = FALSE
convert <- TRUE
skullstrip <- FALSE
plotss = TRUE
regantry <- FALSE
untgantry <- FALSE
runall <- TRUE
useRdcmsort= TRUE
useRdcm2nii= FALSE
removeDups = TRUE
isSorted = NULL
if (ROIformat) isSorted = FALSE
dcm2niicmd = "dcm2nii_2009"
rescale = FALSE
### initial setup
# iid <- length(ids)

iid <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if (is.na(iid)) iid <- 1

id <- ids[iid]
setup(id, study = study)

zeroed <- dir(path=basedir, pattern= ".*Zeroed.*\\.nii\\.gz", 
  recursive=TRUE, full.names=TRUE)
for (ifile in seq_along(zeroed)) {
  system(sprintf('rm "%s"', zeroed[ifile]))
}

# zeroed <- dir(path=homedir, pattern= ":.*.gz", 
# recursive=TRUE, full.names=TRUE)
# for (ifile in seq_along(zeroed)) system(sprintf('rm "%s"', 
# zeroed[ifile]))

# zeroed <- dir(path=homedir, pattern= ":.*.txt", recursive=TRUE, 
# full.names=TRUE)
# for (ifile in seq_along(zeroed)) system(sprintf('rm "%s"', 
# zeroed[ifile]))

# zeroed <- dir(path=homedir, pattern= "'.*.tar.gz", 
# recursive=TRUE, full.names=TRUE)
# for (ifile in seq_along(zeroed)) system(sprintf('rm "%s"', 
# zeroed[ifile]))


#### loop through IDS and convert them to nii, gantry tilted
### 301-520 needs to use Study Date/Time instead of Series Date/Time
# for (iid in 1:length(ids)){

setup(id, study=study)

# for (iid in seq_along(ids)){

  id <- ids[iid]
  print(id)
  setup(id, study = study)
  # source(file.path(progdir, "file_functions.R"))
  dcmsortopt <- ifelse(id =="301-520", '-s ', "")




  if (convert) {
    ### convert the dicoms
  infofile <- file.path(basedir, "Dropout_Information.Rda")
  if (file.exists(infofile)) file.remove(infofile)
  
### time for convertsion
  contime <- NULL
  gf = getfiles(basedir)
# t = iconv("UTF-8","UTF-8//IGNORE",$t);
  
    if (length(gf$files) > 0 | untar){


      contime <- system.time(convert_DICOM(basedir, 
                              verbose=verbose, untar=untar, 
                              useRdcmsort= useRdcmsort, 
                              useRdcm2nii= useRdcm2nii,
                              id = id, 
                              isSorted = isSorted,
                              removeDups=removeDups,
                              dcmsortopt=dcmsortopt, 
                              ROIformat = ROIformat,
                              dcm2niicmd=dcm2niicmd,
                              rescale = rescale))

      # contime <- system.time(convert_DICOM(basedir, progdir, 
      #                         verbose=verbose, untar=untar, 
      #                         useRdcmsort= TRUE, 
      #                         useRdcm2nii= FALSE,
      #                         id = id, 
      #                         dcmsortopt=dcmsortopt, 
      #                         dcm2niicmd=dcm2niicmd))    
      print(contime)

      ## dropout the niis that are not needed
      # lis <- includeMatrix(basedir, dropstring="ungantry", error=TRUE)
      # outs <- lis$outs
      # mis <- lis$mis

      # dropniis <- outs$fname[outs$Takeout]
      # dropniis <- getBase(basename(dropniis), 1)

      # if (length(dropniis) > 0){
      #   dropniis <- file.path(basedir, paste0(dropniis, ".nii.gz"))
      #   for (ifile in dropniis) system(sprintf('rm "%s"', ifile))
      # }

      # save(outs, mis, file = infofile)
    # }
    }
  } ## end of if convert

# }
