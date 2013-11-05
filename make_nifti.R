################################
# Written 2013Oct29
# Author: John Muschelli
# Purpose: Make NIfTI files from DICOMS
# Output: NIfTI images for images and ROIs, and slicethickness image
# Slice thickness image can be used for volume weighting
################################

rm(list=ls())
library(oro.dicom)
basedir <- "/Volumes/DATA_LOCAL/Image_Processing/Test_5"
if (Sys.info()[["user"]] %in% "jmuschel") basedir <- "/dexter/disk2/smart/stroke_ct/ident/Test_5"
dirs <- list.dirs(basedir, recursive=FALSE, full.names=FALSE)
ptpat <- "\\d\\d\\d-(\\d|)\\d\\d\\d"
ids <- grep(paste0(ptpat, "$"), dirs, value=TRUE)


get.tilt <- function(hdr){
  tilt <- as.numeric(hdr$value[ hdr$name %in% "GantryDetectorTilt"])
}


get.slice <- function(hdr){
  tilt <- as.numeric(hdr$value[ hdr$name %in% "SliceThickness"])
}

iid <- 1

iid <- which(grepl("101-308", ids))

# for (iid in 1:length(ids)){
	# testdir <- file.path(basedir, "100-362")
	testdir <- ids[iid]

	setwd(testdir)
	dirs <- list.dirs(path=testdir, recursive=FALSE)
	idir <- 1


	for (idir in 1:2) {
		folname <- basename(dirs[idir])
		print(folname)
		# for (idir in 1:length(dirs)){
		# folname <- "205517_ROI_20110214_1315"
		setwd(file.path(testdir, folname))
		dcm <- readDICOM(path="./", recursive=FALSE, pixelData=FALSE)

		tilts <- sapply(dcm$hdr, get.tilt)
		stopifnot(all(tilts ==0))

		thick <- sapply(dcm$hdr, get.slice)
		if (!all(thick == thick[1])) {
			print(paste0(folname, ": Thickness Wrong"))
		}

		dcm <- readDICOM(path="./", recursive=FALSE)
		# hdr <- dcm$hdr[[1]]

		nim <- dicom2nifti(dcm, rescale=TRUE, reslice=FALSE)
		writeNIfTI(nim, file.path(basedir, folname) )

		### Create Slice Thickness Image
		if (idir %in% 1){
			for (idcm in 1:length(dcm$hdr)){
				thickness <- thick[idcm]
				dims <- dim(dcm$img[[idcm]])
				img <- array(thickness, dim=dims)
				dcm$img[[idcm]] <- img
			}
			nim <- dicom2nifti(dcm, rescale=TRUE, reslice=FALSE)
			filename <- file.path(basedir, 
				"Slice_Thickness",
				folname)
			nim@scl_inter <- 0
			writeNIfTI(nim, filename )				

		}
	}

# }