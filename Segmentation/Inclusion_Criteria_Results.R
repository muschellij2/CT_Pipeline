###################################################################
## This code is for prediction of Image Segmentation of CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
###################################################################
###################################################################
rm(list=ls())
library(plyr)
library(cttools)
library(fslr)
library(ROCR)
library(matrixStats)
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

correct = "none"
options = c("none", "N3", "N4", "N3_SS", "N4_SS",
        "SyN", "SyN_sinc", "Rigid", "Affine")
icorr <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(icorr)) icorr = 1
correct = options[icorr]

keep.obj = ls()


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
    "Affine" = "_Affine")


#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )

outfile = file.path(outdir, "111_Filenames.Rda")
load(file = outfile)


##############################
# Keeping files where predictors exist
##############################
outfiles = nii.stub(basename(fdf$img))
outfiles = paste0(outfiles, "_predictors", adder, ".Rda")
outfiles = file.path(fdf$outdir, outfiles)
stopifnot(file.exists(outfiles))

# load(file = file.path(outdir, "Segmentation_Models.Rda"))
##############################
# Run lmod number of models - not all the models - leave out
##############################
# mod.filename = file.path(outdir, 
# 	paste0("Collapsed_Models", adder, ".Rda"))
# load(mod.filename)

get.pred = 2

for (get.pred in seq(nrow(fdf))){

	iddir = fdf$iddir[get.pred]
	id.outdir = fdf$outdir[get.pred]
	predname = nii.stub(basename(fdf$img[get.pred]))
	predname = file.path(id.outdir, 
		paste0(predname, "_predictors", adder, ".Rda"))
	load(predname)
	df = img.pred$df
	nim = img.pred$nim
	keep.ind = img.pred$keep.ind
	rm(list="img.pred")
    for (i in 1:3) gc()	
	df$include = df$value >= 30 & df$value <= 100


    fname = file.path(outdir, 
        paste0("Aggregate_data_cutoffs", adder, ".Rda"))

    load(file = fname)
    keepnames = colnames(est.cutoffs)
    include = rep(TRUE, length=nrow(df))
    for (icut in keepnames){
    	qcuts = est.cutoffs[, icut]
    	include = include & 
    		(df[, icut] >= qcuts[1] & df[, icut] <= qcuts[2])
    }

	df$include.all = include

		# df$dist_centroid <= 75

	df$in0100 = df$value >= 0 & df$value <= 100
	# df$in20_85 = df$value >= 20 & df$value <= 85
	df$mask = df$mask > 0
	
	######################################
	# Get volume of ROI
	######################################	
	not0100 = sum(df$Y[ !df$in0100 ])

	######################################
	# Keep all ROI = 1, even if not inmask
	######################################	
	roi.not.in = which(df$Y == 1)
	roi.not.in = roi.not.in[!(roi.not.in %in% keep.ind)]
	keep.ind = sort(c(keep.ind, roi.not.in))


	#### need this because the length of df has changed
	roi.not.in = which(keep.ind %in% roi.not.in)
	nroi.not.in = length(roi.not.in)

	df = df[keep.ind,]

	pct.out = sum(df$Y[! df$include.all]) / sum(df$Y)
	pct.reduced = sum(1-df$Y[!df$include.all]) / sum(1-df$Y)

	pct.30out = sum(df$Y[! df$include]) / sum(df$Y)
	pct.30reduced = sum(1-df$Y[!df$include]) / sum(1-df$Y)

	print(pct.out)
	print(pct.30out)
	print(pct.reduced)
	print(pct.30reduced)
	print(get.pred)
	predname = nii.stub(basename(fdf$img[get.pred]))
	predname = file.path(id.outdir, 
		paste0(predname, "_cutoff_results", adder, ".Rda"))	
	save(pct.out, pct.reduced, nroi.not.in,
		pct.30out, 
		pct.30reduced,
		file = predname)
} 