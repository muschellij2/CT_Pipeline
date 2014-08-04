#################################
# Creates map of top p-values from unadusted model
# Author: John Muschelli
#################################
rm(list=ls())
library(oro.nifti)
library(plyr)
library(scales)
library(RColorBrewer)
library(data.table)
library(cttools)
library(fslr)
library(smallPDF)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
# source(file.path(progdir, "convert_DICOM.R"))
# source(file.path(progdir, "fslhd.R"))
basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")


template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)


whichdir = "reoriented"
outcome = "NIHSS"
adder = paste0(outcome, "_")
if (outcome == "NIHSS"){
	adder = ""
}

ID = "100-362"

rerun = FALSE
dtemp = dim(temp)
	### go to id directory
	# if (verbose) print(ID)
	iddir = file.path(basedir, ID)
	resdir = file.path(iddir, whichdir)

	imgs = file.path(resdir, 
		"bws100-362_20100126_1926_CT_2_CT_ROUTINEROI.nii")	
	rawimgs = file.path(resdir, 
		"w100-362_20100126_1926_CT_2_CT_ROUTINE.nii")
	ssimgs = file.path(iddir, "Skull_Stripped",
		"100-362_20100126_1926_CT_2_CT_ROUTINE_SS_0.01.nii.gz")
	baseimgs = file.path(iddir, 
		"100-362_20100126_1926_CT_2_CT_ROUTINE.nii.gz")
	temprois = file.path(resdir,  
		"bws100-362_20100126_1926_CT_2_CT_ROUTINEROI.nii")
	rawrois = file.path(rootdir, "ROI_data", ID,
		"100-362_20100126_1926_CT_2_CT_ROUTINEROI.nii.gz")	

	id = basename(rawimgs)
	id = gsub("\\.gz$", "", id)
	iid = id = gsub("\\.nii$", "", id)

	id = gsub("^(w|affine(9|12)_)", "", iid)
	rawpng= file.path(resdir, paste0("raw_", id, ".png"))
	opng = file.path(resdir, paste0("roi_", id, ".png"))
	sspng = file.path(resdir, paste0("ss_", id, ".png"))
	rawroipng = file.path(resdir, paste0("native_", id, ".png"))

	verbose= TRUE

	# if (rerun){


		rawimg = readNIfTI(rawimgs)
		rawimg[is.nan(rawimg) | is.na(rawimg)] = 0
		# rawimg[ rawimg > 100 | rawimg < 0] = 0
		rawmask = rawimg >= 0 & rawimg <= 100


		img = readNIfTI(imgs)
		img[is.nan(img) | is.na(img)] = 0
		tt = temp
		tt[ tt >= 99 | tt < 0] = 0
		tt[ img > 0 ] = 100

		## get cross hair information
		chair = round(colMeans(which(img > 0, arr.ind=TRUE)))

		chair[is.nan(chair)] = dtemp[is.nan(chair)]/2
		png(opng, type="cairo", res=600, 
			height=7,width=7, units="in")

			orthographic(tt, 
				col= c(gray(0), gray(1:60/60), hotmetal(1)), 
				xyz=chair, text="D", text.cex = 8)

		dev.off()

		png(rawpng, type="cairo", res=600, 
			height=7,width=7, units="in")

			# orthographic(rawimg,  
			# 	xyz=chair)
			mask.overlay(rawimg, img, xyz=chair, 
				col.y=alpha("red", 0.25), 
				text="C", text.cex = 8)

		dev.off()


	#### Raw data

		ssimg = readNIfTI(ssimgs)
		ssimg[is.nan(ssimg) | is.na(ssimg)] = 0
		dimg = dim(ssimg)

		img = readNIfTI(rawrois)
		img[is.nan(img) | is.na(img)] = 0

		chair = round(colMeans(which(img > 0, arr.ind=TRUE)))

		chair[is.nan(chair)] = dimg[is.nan(chair)]/2
		
		png(sspng, type="cairo", res=600, 
			height=7,width=7, units="in")

			mask.overlay(ssimg, img, xyz=chair, 
				col.y=alpha("red", 0.25), 
				text="B", text.cex = 8)

		dev.off()

		baseimg = readNIfTI(baseimgs)
		baseimg[is.nan(baseimg) | is.na(baseimg)] = 0

		png(rawroipng, type="cairo", res=600, 
			height=7,width=7, units="in")

			mask.overlay(baseimg, img, xyz=chair, 
				col.y=alpha("red", 0.25), 
				text="A", text.cex = 8)

		dev.off()
		
		cat(paste0("rsync --progress $jenig:", 
			c(rawroipng, sspng, rawpng, opng),
			" ~/CT_Registration/CT_Pipeline/"), sep="\n")




