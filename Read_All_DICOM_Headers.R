################################
# Written 2013Oct29
# Author: John Muschelli
# Purpose: Make NIfTI files from DICOMS
# Output: NIfTI images for images and ROIs, and slicethickness image
# Slice thickness image can be used for volume weighting
################################

rm(list=ls())
library(oro.dicom)
library(oro.nifti)
library(plyr)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"  
}
rootprogdir = file.path(rootdir, "programs")
## need to make package
source(file.path(rootprogdir, "fslhd.R"))
progdir = file.path(rootprogdir, "Test_Registration")
basedir = file.path(rootdir, "MISTIE")
dirs <- list.dirs(basedir, recursive=FALSE, full.names=FALSE)
ptpat <- "\\d\\d\\d-(\\d|)\\d\\d\\d"
ids <- grep(paste0(ptpat, "$"), dirs, value=TRUE)

iid =1

all.hdrs = list()
for (iid in seq_along(ids)){
	testdir <- file.path(basedir, ids[iid])

	setwd(testdir)
	### get all dcm files
	dcms <- list.files(path=testdir, pattern= ".dcm$",
		recursive=TRUE, full.names=TRUE)

	print(ids[iid])
	# ifile = 1;
	### read in EVERY HEADER from this 
	hdrl = llply(dcms, 
		function(x) {
			rereadDICOMFile(x, pixelData=FALSE)$hdr
		}, .progress = "text")
	
	names(hdrl) = dcms

	dcmtables = dicomTable(hdrl)
	save(dcmtables, file=file.path(testdir, "All_Headers.Rda"))
}