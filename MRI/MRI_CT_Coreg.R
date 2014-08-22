rm(list=ls())
library(cttools)
library(fslr)
library(plyr)
options(matlab.path='/Applications/MATLAB_R2013b.app/bin')

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


setup(id, study=study)

# for (iid in seq_along(ids)){

id <- ids[iid]
print(id)
setup(id, study = study)
  
ctfile = file.path(basedir, 
  "100-365_20100217_2338_CT_2_CT_Head_Spiral_Thresh.nii.gz")
# mrifile = file.path(basedir, 
#   "100-365_20100218_1837_MR_3_BRAIN_T2.nii.gz")
mrifile = file.path(basedir, 
  "100-365_20100218_1829_MR_2_BRAIN_T1.nii.gz")

q = 0.5
ct_thresh = -100
opts = "-v"
dof = 6
outfile = paste0(nii.stub(mrifile), "_Coreg")
omat = paste0(outfile, ".mat")
retimg = FALSE
modality = "T1"

mri_thresh = function(mrifile= mrifile, q=q, modality){
  if (modality %in% c("T2", "T1")){
    quant = quantile(c(mrifile), probs = q, na.rm=TRUE)
    mrifile[which(mrifile < quant)] = NA
    mrifile = cal_img(mrifile)
    return(mrifile)
  }
}

# MR_CT_flirt = function(ctfile, 
# mrifile, q = 0.5, 
# ct_thresh = -100, dof = 6, outdir = ...){
  
  mrifile = check_nifti(mrifile, reorient=FALSE)
  ctfile = check_nifti(ctfile, reorient=FALSE)

  q =  q[1]
  mrifile = mri_thresh(mrifile= mrifile, q=q, modality=modality)

  ctfile[ctfile < ct_thresh ] = NA


  res = flirt(infile = mrifile, 
    reffile = ctfile, 
    omat = omat,
    outfile = outfile,
    retimg = retimg,
    dof = dof, opts = opts, ...)

# }