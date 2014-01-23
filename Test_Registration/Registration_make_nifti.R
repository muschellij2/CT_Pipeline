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
basedir = file.path(rootdir, "Test_Registration")
dirs <- list.dirs(basedir, recursive=FALSE, full.names=FALSE)
ptpat <- "\\d\\d\\d-(\\d|)\\d\\d\\d"
ids <- grep(paste0(ptpat, "$"), dirs, value=TRUE)

get.nifti.header = function(nim){
	sn = slotNames(nim)
	sn = sn[ sn != ".Data"]
	hdr = sapply(sn, slot, object=nim)
}

# Sys.setenv("EXIT_CODE"="100")
# Sys.setenv("EXIT_STATUS"="100")


get.tilt <- function(hdr){
  tilt <- as.numeric(hdr$value[ hdr$name %in% "GantryDetectorTilt"])
}


get.slice <- function(hdr){
  tilt <- as.numeric(hdr$value[ hdr$name %in% "SliceThickness"])
}

iid <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if (is.na(iid)) iid <- 11

# ids <- ids[iid]

# iid <- which(grepl("101-308", ids))

# for (iid in 1:length(ids)){
	# testdir <- file.path(basedir, "100-362")
	testdir <- file.path(basedir, ids[iid])

	setwd(testdir)
	dirs <- list.dirs(path=testdir, recursive=FALSE, full.names=TRUE)
	idir <- 1


	for (idir in 1:2) {
		folname <- basename(dirs[idir])
		print(folname)
		# for (idir in 1:length(dirs)){
		# folname <- "205517_ROI_20110214_1315"
		setwd(dirs[idir])
		dcm <- readDICOM(path="./", 
			recursive=FALSE, 
			pixelData=FALSE)

		tilts <- sapply(dcm$hdr, get.tilt)
		stopifnot(all(tilts ==0))

		thick <- sapply(dcm$hdr, get.slice)
		if (!all(thick == thick[1])) {
			print(paste0(folname, ": Thickness Wrong"))
		}

		dcm <- readDICOM(path="./", recursive=FALSE)
		# hdr <- dcm$hdr[[1]]

		dcmtable = dicomTable(dcm$hdr)
		keepcols = grepl("RescaleIntercept|RescaleSlope", 
			colnames(dcmtable))
		dcmtab = dcmtable[, keepcols]
		stopifnot(ncol(dcmtab) == 2)
		colnames(dcmtab) = c("intercept", "slope")


		iimg = 1;
		### need to manually change values - windowing problems FOV
		for (iimg in 1:length(dcm$img)){
			inter = as.numeric(dcmtab$intercept[iimg])
			slope = as.numeric(dcmtab$slope[iimg])
			dcm$img[[iimg]] = (dcm$img[[iimg]] + inter)/slope
			x = dcm$img[[iimg]]
			print(range(dcm$img[[iimg]]))
			dcm$img[[iimg]][x < -1024] = -1024
			dcm$img[[iimg]][x > 3000] = 3000
		}
		nim <- dicom2nifti(dcm, rescale=FALSE, reslice=FALSE)
		nim@scl_slope = 1
		if (idir == 2){
			nim <- nim + nim@scl_inter
			nim@cal_max <- nim@cal_max + nim@scl_inter
			nim@cal_min <- nim@cal_min + nim@scl_inter
			nim@scl_inter <- 0
		}
		full.name = file.path(basedir, "RawNIfTI", folname)
		writeNIfTI(nim, 
			full.name, 
			gzipped=FALSE )
		nifti.header = get.nifti.header(nim)
		fsl.header = fslhd.parse(fslhd(full.name))
		save(dcmtab, nifti.header, fsl.header,
			file=paste0(full.name, "_Header.Rda"))		
		## copy over for reorientation
		reo = file.path(basedir, "reoriented", folname)
		stopifnot(file.copy(paste0(full.name, '.nii'), 
			paste0(reo, '.nii'), overwrite=TRUE))


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
			writeNIfTI(nim, filename, gzipped=TRUE)				

		}
	}

# }

Sys.setenv("exit_code"="0")
# Sys.setenv("exit_code"="0")
