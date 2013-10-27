rm(list=ls())
library(oro.dicom)
basedir <- "/Volumes/Seagate Backup Plus Drive/Image_Processing/Test_5"
testdir <- file.path(basedir, "301-520")
setwd(testdir)
dirs <- list.dirs(path=testdir, recursive=FALSE)
idir <- 2
folname <- basename(dirs[idir])
# for (idir in 1:length(dirs)){
# folname <- "205517_ROI_20110214_1315"
setwd(file.path(testdir, folname))
dcm <- readDICOM(path="./", recursive=FALSE)
# hdr <- dcm$hdr[[1]]
get.tilt <- function(hdr){
  tilt <- as.numeric(hdr$value[ hdr$name %in% "GantryDetectorTilt"])
}
tilts <- sapply(dcm$hdr, get.tilt)
stopifnot(all(tilts ==0))
nim <- dicom2nifti(dcm, rescale=TRUE, reslice=FALSE)
writeNIfTI(nim, file.path(testdir, folname) )