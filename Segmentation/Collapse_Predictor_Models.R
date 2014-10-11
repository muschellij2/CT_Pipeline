###################################################################
## This code is for prediciton of Image Segmentation of CT
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

correct = "N3"
options = c("none", "N3", "N4", "N3_SS", "N4_SS",
        "SyN", "SyN_sinc", "Rigid", "Affine")

for (correct in options){
    rm(list="all.df")
    for (i in 1:3) gc()
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
	non.aggmods = lmod = 10
	fdf.run = fdf[seq(lmod), ]
	#### adding aggregate model
	lmod = lmod +1
	imod = 6
	runpreds = 1:nrow(fdf)
	res = matrix(NA, nrow = lmod, ncol = nrow(fdf))
	rownames(res) = paste0("mod", seq(lmod))
	# colnames(res) = paste0("pred", seq(nrow(fdf)))
	### number of iterations of predictions from models
	nextra = 6
	sres = res
	vol.data = matrix(NA, nrow= length(runpreds), ncol =lmod+1+nextra)
	vol.sdata = vol.data
	res = matrix(NA, nrow= length(runpreds), ncol =lmod+nextra)
	sres = res
	# get.pred = 110
	colnames(vol.data) = colnames(vol.sdata) = 
		c("truth", paste0("model", seq(lmod)), 
			"mean", "median", "max", "min", "prod", "gmean")
	colnames(res) = colnames(sres) = colnames(vol.data)[-1]
	colnames(vol.data)[lmod] = colnames(vol.sdata)[lmod] =
		colnames(res)[lmod] = "mod_agg"

	fdf.run = rbind(fdf.run, NA)
	fdf.run$img = as.character(fdf.run$img)
	fdf.run$img[lmod] = "Aggregate"
	fdf.run$outdir[lmod] = outdir

	all.mods = llply(seq(lmod), function(imod){
			mod.outdir = fdf.run$outdir[imod]
			moddname = nii.stub(basename(fdf.run$img[imod]))
			moddname = file.path(mod.outdir, 
				paste0(moddname, "_models", adder, ".Rda"))

			load(moddname)
			mod = mods$mod
			return(mod)
	}, .progress = "text")
	filename = file.path(outdir, 
		paste0("Collapsed_Models", adder, ".Rda"))
	save(all.mods, res, vol.data, vol.sdata, sres, fdf.run,
		runpreds,
		lmod, non.aggmods, file=filename)
	print(correct)
}