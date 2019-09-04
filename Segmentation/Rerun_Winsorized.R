###########################################################
## This code is for Creating Predictors for each individual
##
## Author: John Muschelli
## Last updated: May 20, 2014
############################################################
############################################################
rm(list=ls())
library(plyr)
library(cttools)
library(devtools)
library(ROCR)
library(fslr)
library(mgcv)
library(getopt)
library(psych)
library(methods)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")
outdir = file.path(basedir, "results")

template.file = file.path(tempdir, "scct_unsmooth.nii.gz")
ss.tempfile = file.path(tempdir, "Skull_Stripped",
    "scct_unsmooth_SS_First_Pass_0.1.nii.gz")


mean.file = file.path(tempdir, "Mean_Image.nii.gz")
sd.file = file.path(tempdir, "SD_Image.nii.gz")

# mean.img = readNIfTI(mean.file)
# sd.img = readNIfTI(sd.file)


segdir = file.path(progdir, "Segmentation")
source(file.path(segdir, "performance_functions.R"))

spec = matrix(c(
	'correct', 'c', 1, "character",
	'rerun'   , 'r', 1, "logical",
	'overwrite'  , 'o', 1, "logical",
	'interpolator'  , 'i', 1, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

correct = opt$correct
rerun = opt$rerun
overwrite = opt$overwrite
interpolator = opt$interpolator

if (is.null(rerun)){
	rerun = TRUE
}
if (is.null(overwrite)){
	overwrite = TRUE
}
if (is.null(interpolator)){
	interpolator = "Linear"
}
print(opt)
# args<-commandArgs(trailingOnly = TRUE)
# print(args)
# correct = args[1]
# rerun = as.logical(args[2])
# overwrite = as.logical(args[3])
# print(rerun)
# print(overwrite)

# correct = "Rigid_sinc"
# options = c("none", "N3", "N4", "N3_SS", "N4_SS",
#         "SyN", "SyN_sinc", "Rigid", 
# "Affine", "Rigid_sinc", 
#         "Affine_sinc")
# options = c("none", "N3", "N4", "N3_SS", "N4_SS", 
# 		"Rigid", "Rigid_sinc")
options = c("none", "N3_SS", "N4_SS", 
		"Rigid", "Rigid_sinc")
# options = c("Rigid", "Rigid_sinc")

#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )

outfile = file.path(outdir, "111_Filenames.Rda")
load(file = outfile)


# iscen <- as.numeric(Sys.getenv("SGE_TASK_ID"))
# if (is.na(iscen)) iscen = 1


# scenarios = expand.grid(iimg = seq(nrow(fdf)), 
# 	correct = options, stringsAsFactors = FALSE )
# iimg = scenarios$iimg[iscen]
# correct = scenarios$correct[iscen]

iimg <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(iimg)) iimg = 1

runx = x = fdf[iimg,]

# for (correct in options){
    
    correct = match.arg(correct, options)

	adder = switch(correct, 
		"none"= "",
		"N3"="_N3",
		"N4" = "_N4",
		"N3_SS" = "_N3_SS",
		"N4_SS" = "_N4_SS", 
		"SyN" = "_SyN",
		"SyN_sinc" = "_SyN_sinc",
		"Rigid" = "_Rigid",
		"Affine" = "_Affine",
		"Rigid_sinc" = "_Rigid_sinc",
		"Affine_sinc" = "_Affine_sinc")

	# fdf = fdf[1:10,]
# run_model = function(x, fpr.stop = .1){
	fname= nii.stub(basename(x$img))
	fname = paste0(fname, "_predictors", adder, ".Rda")
	outfile = file.path(x$outdir, fname)
	fname = switch(correct,
		"none" = x$img,
		"N3" = x$n3img,
		"N4" = x$n4img,
		"N3_SS" = x$n3ssimg,
		"N4_SS" = x$n4ssimg,
		"SyN" = x$synssimg,
		"SyN_sinc" = x$sinc_synssimg,
		"Rigid" = x$rig_ssimg,
		"Affine" = x$aff_ssimg,
		"Rigid_sinc" = x$sinc_rig_ssimg,
		"Affine_sinc" = x$sinc_aff_ssimg		
		)
	mask.fname = switch(correct,
		"none" = x$mask,
		"N3" = x$mask,
		"N4" = x$mask,
		"N3_SS" = x$mask,
		"N4_SS" = x$mask,
		"SyN" = x$synssmask,
		"SyN_sinc" = x$sinc_synssmask,
		"Rigid" = x$rig_ssmask,
		"Affine" = x$aff_ssmask,
		"Rigid_sinc" = x$sinc_rig_ssmask,
		"Affine_sinc" = x$sinc_aff_ssmask
		)
	roi.fname = switch(correct,
		"none" = x$roi,
		"N3" = x$roi,
		"N4" = x$roi,
		"N3_SS" = x$roi,
		"N4_SS" = x$roi,
		"SyN" = x$synssroi,
		"SyN_sinc" = x$sinc_synssroi,
		"Rigid" = x$rig_ssroi,
		"Affine" = x$aff_ssroi,
		"Rigid_sinc" = x$sinc_rig_ssroi,
		"Affine_sinc" = x$sinc_aff_ssroi
		)	
	print(correct)


	xxx = load(file=outfile)
	cn = colnames(img.pred$df)
	if (!("win_z" %in% cn)){
		ofile = x$img
		if (correct %in% c("Rigid", "Rigid_sinc")){
			ofile = x$ss
		}

		ofile = file.path( x$outdir, 
			paste0(nii.stub(ofile, bn=TRUE), 
				adder, "_win_z", ".nii.gz"))
		img = readNIfTI(ofile, reorient = FALSE)
		img.pred$df$win_z = c(img)
		save(img.pred, file=outfile, compress = TRUE)
	}
