###################################################################
## This code is for prediciton of Image Segmentation of CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
###################################################################
###################################################################
rm(list=ls())
library(ggplot2)
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

correct = "SyN"
options = c("none", "N3", "N4", "N3_SS", "N4_SS",
		"SyN", "SyN_sinc", "Rigid", "Affine")
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
		"Affine" = "_Affine")


	#### load voxel data
	outfile = file.path(outdir, "Voxel_Info.Rda")
	load(file=outfile )

	outfile = file.path(outdir, "111_Filenames.Rda")
	load(file = outfile)


	makedir = sapply( fdf$outdir, function(x) {
		if (!file.exists(x)){
			dir.create(x, showWarnings =FALSE)
		}
	})
	irow = 1
	x = fdf[irow,]

	##############################
	# Keeping files where predictors exist
	##############################
	outfiles = nii.stub(basename(fdf$img))
	outfiles = paste0(outfiles, "_predictors", adder, ".Rda")
	outfiles = file.path(fdf$outdir, outfiles)
	stopifnot(file.exists(outfiles))

	get.pred <- as.numeric(Sys.getenv("SGE_TASK_ID"))
	if (is.na(get.pred)) get.pred = 16
	x = fdf[get.pred,]


	pdfname = file.path(outdir, 
	    paste0("ROI_Histogram", adder, ".pdf"))
	pdf(pdfname)
	for (get.pred in seq(nrow(fdf))){

		iddir = fdf$iddir[get.pred]
		outdir = fdf$outdir[get.pred]
		img.stub = nii.stub(fdf$img[get.pred], bn=TRUE)
		predname = nii.stub(basename(fdf$img[get.pred]))
		predname = file.path(outdir, 
			paste0(predname, "_predictors", adder, ".Rda"))
		load(predname)
		df = img.pred$df
		nim = img.pred$nim
		keep.ind = img.pred$keep.ind
		df$in0100 = df$value >= 0 & df$value <= 100
		df$mask = df$mask > 0


		######################################
		# Keep all ROI = 1, even if not inmask
		######################################	
		roi.not.in = which(df$Y == 1)
		df = df[roi.not.in,]

		hist(df$value, main=img.stub)
		print(get.pred)
	}
	dev.off()

}