rm(list=ls())
library(cttools)
library(fslr)
library(plyr)
options(matlab.path='/Applications/MATLAB_R2013b.app/bin')

rootdir = "~/Lung/data"
rootdir = path.expand(rootdir)

iddirs = list.dirs(path=rootdir, full.names=TRUE, recursive=FALSE)
ids = basename(iddirs)

iid = 1
for (iid in seq_along(ids)){
    basedir = iddirs[iid]

    #### setting up if things are on the cluster or not
    ROIformat = FALSE


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

    ### time for convertsion
      contime <- NULL
      gf = getfiles(basedir)
      


    contime <- system.time(convert_DICOM(basedir, 
                            verbose=verbose, untar=untar, 
                            useRdcmsort= useRdcmsort, 
                            useRdcm2nii= useRdcm2nii,
                            id = id, 
                            isSorted = isSorted,
                            removeDups=removeDups,
                            dcmsortopt=dcmsortopt, 
                            ROIformat = ROIformat,
                            dcm2niicmd=dcm2niicmd))

}