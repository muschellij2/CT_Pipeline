##########################################################
## This code is for flipping images and taking differences
##
## Author: John Muschelli
## Last updated: May 20, 2014
###########################################################
###########################################################
rm(list=ls())
library(plyr)
library(cttools)
library(ROCR)
library(matrixStats)
library(reshape2)
library(ggplot2)
library(fslr)
library(extrantsr)
library(ANTsR)
library(getopt)
# library(car)
library(GGally)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
segdir = file.path(progdir, "Segmentation")
basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

template.file = file.path(tempdir, 
	"scct_unsmooth_SS_0.01.nii.gz")
segdir = file.path(progdir, "Segmentation")
source(file.path(segdir, "performance_functions.R"))

outdir = file.path(basedir, "results")
correct = "none"

# options = c("none", "N3_SS", "N4_SS", 
# 		"Rigid", "Rigid_sinc")
options = c("none", "Rigid")

# types = c("", "_include", "_zval", "_zval2")
# types = c("_zval2", "_zval_all", '_zval2_medztemp')
types = "_zval2"
# types = "_zval2"
# "_include_all", 
# types = "_include_all"
type = types[1]

spec = matrix(c(
	'correct', 'c', 1, "character",
	'rerun'   , 'r', 1, "logical"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

correct = opt$correct
rerun = opt$rerun

if (is.null(rerun)){
	rerun = TRUE
}
# for (correct in options){

print(correct)
correct = match.arg(correct, options)

#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )

outfile = file.path(outdir, 
	"111_Filenames_with_volumes_stats.Rda")
load(file = outfile)

##########################
# Select subject
##########################
iimg <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(iimg)) iimg = 17

runx = x = fdf[iimg,]

##########################
# Make fnames
##########################
fname = switch(correct,
	"none" = x$ssimg,
	"Rigid" = x$rig_ssimg	
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


##########################
# Prep outfiles
##########################
outfile = file.path(x$outdir, 
    paste0(nii.stub(fname, bn=TRUE), 
	"_Flipped_Difference.nii.gz"))

if (!all(file.exists(outfile)) | rerun) {

	mask = fslbin(
		fslthresh(mask.fname, thresh = 0.5)
	)

	verbose = TRUE
	mask.outfile = NULL
	typeofTransform = "Rigid"
	t1.outfile = NULL
	img = readnii(fname)

    flipper = reg_flip( img, 
        mask = mask,
        template.file = template.file, 
        typeofTransform = "Affine", 
        interpolator = "LanczosWindowedSinc")
    fres = flipper$t1

	template = readnii(template.file)
	outprefix = tempfile()
	##########################
	# register to template
	##########################
    res = ants_regwrite(
    	filename = fname,
    	outfile = NULL, 
    	retimg = TRUE, 
        template.file = template.file, 
        typeofTransform = "Affine", 
		interpolator = "LanczosWindowedSinc", 
		outprefix = outprefix, 
        remove.warp = TRUE, 
        verbose = TRUE)	
    ##########################
    # Flip image
    ##########################
    fres = fslswapdim(res, 
    	retimg = TRUE, a = "-x")
    ##########################
    # Take difference
    diffimg = res - fres

	##########################
	# revert space
	##########################
    inv.trans = paste0(outprefix, "0GenericAffine.mat")
    fixed = oro2ants(fname)
    moving = oro2ants(diffimg)
    diffimg.native = antsApplyTransforms(
    	fixed = fixed, 
    	moving = moving, 
      transformlist = inv.trans, 
      interpolator = "LanczosWindowedSinc", 
      whichtoinvert = 1)
	##########################
	# Mask flipped difference image
	##########################    
    dimg = mask_img(
    	ants2oro(diffimg.native),
    	mask)

	##########################
	# Write out image
	##########################    
    writeNIfTI(dimg, 
    	filename = nii.stub(outfile))
}
