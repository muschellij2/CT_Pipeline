##############################################
## This code is for Creating Predictors 
# for each individual
##
## Author: John Muschelli
## Last updated: May 20, 2014
##############################################
##############################################
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

# template.file = file.path(tempdir, 
# 	"scct_unsmooth.nii.gz")
# ss.tempfile = file.path(tempdir, "Skull_Stripped",
#     "scct_unsmooth_SS_First_Pass_0.1.nii.gz")


# mean.file = file.path(tempdir, "Mean_Image.nii.gz")
# sd.file = file.path(tempdir, "SD_Image.nii.gz")

# mean.img = readNIfTI(mean.file)
# sd.img = readNIfTI(sd.file)

segdir = file.path(progdir, "Reseg")
source(file.path(segdir, 
	"Reseg_performance_functions.R"))

correct = "Rigid"

spec = matrix(c(
	'rerun'   , 'r', 1, "logical",
	'overwrite'  , 'o', 1, "logical",
	'interpolator'  , 'i', 1, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

rerun = opt$rerun
overwrite = opt$overwrite
interpolator = opt$interpolator

if (is.null(rerun)){
	rerun = FALSE
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
options = c("Rigid")
# options = c("Rigid", "Rigid_sinc")

#### load image information
outfile = file.path(outdir, 
	"Reseg_111_Filenames.Rda")
load(file = outfile)


iimg <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(iimg)) iimg = 21

runx = x = fdf[iimg,]

# for (correct in options){
    
    correct = match.arg(correct, options)

	adder = switch(correct, 
		"none"= "",
		"Rigid" = "_Rigid")

	# fdf = fdf[1:10,]
# run_model = function(x, fpr.stop = .1){
	fname= nii.stub(basename(x$img))
	fname = paste0("Reseg_", fname, 
		"_predictors", adder, 
		".Rda")
	outfile = file.path(x$outdir, 
		fname)
	fname = switch(correct,
		"none" = x$img,
		"Rigid" = x$rig_ssimg	
		)
	mask.fname = switch(correct,
		"none" = x$mask,
		"Rigid" = x$rig_ssmask
		)
	roi.fname = switch(correct,
		"none" = x$roi,
		"Rigid" = x$rig_ssroi
		)	
	print(correct)

	roi = (readnii(roi.fname) > 0.5)*1
	mask = readnii(mask.fname) > 0.5
	if (!file.exists(outfile) | rerun){
	  	system.time({
	  		img.pred = make_predictors(
	  		img= fname, 
	  		stub = nii.stub(fname, bn = TRUE),
	  		mask = mask, 
	  		roi = roi,
	  		save_imgs = TRUE, 
	  		outdir = x$outdir,   
	  		overwrite = overwrite, 
	  		verbose= TRUE)
	  	})
		img.pred$df$Y = as.numeric(
			img.pred$df$Y > 0.5
			)	  	
		save(img.pred, file=outfile, 
			compress = TRUE)
	} 

	img.pred$df$Y = as.numeric(
		img.pred$df$Y > 0.5
		)
	
	# cut = c(-10, 10)
	# img.pred$df$zscore_template[
	# 	img.pred$df$zscore_template >= cut[2]
	# 	] = cut[2]
	# img.pred$df$zscore_template[
	# 	img.pred$df$zscore_template <= cut[1]
	# 	] = cut[1]
	# img.pred$df$Y = (img.pred$df$Y > 0.5)*1
	# save(img.pred, file=outfile, compress = TRUE)

	# stopifnot(all(img.pred$df$Y %in% c(0, 1)))
# }
	