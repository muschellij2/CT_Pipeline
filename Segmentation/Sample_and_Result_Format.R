# 
## This code is for collapsing predictor models
##
## Author: John Muschelli
## Last updated: May 20, 2014
#########################################################
#########################################################
rm(list=ls())
library(plyr)
library(fslr)
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
# options = c("none", "N3", "N4", "N3_SS", "N4_SS",
#         "SyN", "SyN_sinc", "Rigid", "Affine", 
# "Rigid_sinc", 
#         "Affine_sinc")
options = c("none", "N3_SS", "N4_SS", 
      "Rigid", "Rigid_sinc")

icorr <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(icorr)) icorr = 1
correct = options[icorr]


for (correct in options){

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

	#### load voxel data
	outfile = file.path(outdir, "Voxel_Info.Rda")
	load(file=outfile )

	outfile = file.path(outdir, "111_Filenames.Rda")
	load(file = outfile)


	##############################
	# Keeping files where predictors exist
	##############################
	# outfiles = nii.stub(basename(fdf$img))
	# outfiles = paste0(outfiles, "_predictors", 
	# adder, ".Rda")
	# outfiles = file.path(fdf$outdir, outfiles)
	# stopifnot(file.exists(outfiles))

	# load(file = file.path(outdir, 
		# "Segmentation_Models.Rda"))
	##############################
	# Run lmod number of models - 
	# not all the models - leave out
	##############################


	run.ind = which(fdf$group == "Train")
  # run.ind = seq(lmod)
	fdf.run = fdf[run.ind , ]
	non.aggmods = nrow(fdf.run)
	#### adding aggregate model
	lmod = non.aggmods + 1
	imod = 6
	runpreds = 1:nrow(fdf)
	res = matrix(NA, nrow = lmod, ncol = nrow(fdf))
	rownames(res) = paste0("mod", seq(lmod))
	# colnames(res) = paste0("pred", seq(nrow(fdf)))
	### number of iterations of predictions from models
	nextra = 8
	sres = res
	vol.data = matrix(NA, nrow= length(runpreds), 
		ncol = lmod+1+nextra)
	vol.sdata = vol.data
	res = matrix(NA, nrow= length(runpreds), 
		ncol =lmod+nextra)
	sres = res
	# get.pred = 110
	colnames(vol.data) = colnames(vol.sdata) = 
		c("truth", paste0("model", seq(lmod)), 
			"mean", "median", "max", "min", 
			"prod", "gmean", 
			"gam", "rf")
	colnames(res) = colnames(sres) = 
		colnames(vol.data)[-1]
	colnames(sres)[lmod] = colnames(res)[lmod] = 
		"mod_agg"
	colnames(vol.data)[lmod+1] = 
		colnames(vol.sdata)[lmod+1] = 
		"mod_agg"

	fdf.run = rbind(fdf.run, NA)
	fdf.run$img = as.character(fdf.run$img)
	fdf.run$img[lmod] = "Aggregate"
	fdf.run$outdir[lmod] = outdir

	filename = file.path(outdir, 
		paste0("Result_Formats", adder, ".Rda"))
	save(res, vol.data, vol.sdata, sres, fdf.run,
		runpreds, run.ind,
		lmod, non.aggmods, file=filename)
	print(adder)

}

