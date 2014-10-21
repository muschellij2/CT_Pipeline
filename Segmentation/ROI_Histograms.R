###################################################################
## This code is for Histograms of ROI
##
## Author: John Muschelli
## Last updated: May 20, 2014
###################################################################
###################################################################
rm(list=ls())
library(ggplot2)
library(fslr)
library(plyr)
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

# correct = "SyN"
options = c("none", "N3", "N4", "N3_SS", "N4_SS",
		"SyN", "SyN_sinc", "Rigid", "Affine")
icorr <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(icorr)) icorr = 1
correct = options[icorr]

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

	get.pred = 16
	x = fdf[get.pred,]


	pdfname = file.path(outdir, 
	    paste0("ROI_Histogram", adder, ".pdf"))
	pdf(pdfname)
	for (get.pred in seq(nrow(fdf))){

		iddir = fdf$iddir[get.pred]
		id.outdir = fdf$outdir[get.pred]
		img.stub = nii.stub(fdf$img[get.pred], bn=TRUE)
		predname = img.stub
		predname = file.path(id.outdir, 
			paste0(predname, "_predictors", adder, ".Rda"))
		load(predname)
		df = img.pred$df
		nim = img.pred$nim
		keep.ind = img.pred$keep.ind
		# df$in0100 = df$value >= 0 & df$value <= 100
		df$mask = df$mask > 0

		######################################
		# Keep all ROI = 1, even if not inmask
		######################################			
		roi.in = which(df$Y == 1)
		roi.not.in = which(df$Y != 1)

		roi.not.inc = roi.in[!(roi.in %in% keep.ind)]
		keep.ind = sort(c(keep.ind, roi.not.inc))

		##### need to subset relevant data
		df = df[keep.ind, ]

		######################################
		# Get indices of correct data
		######################################	
		roi.in = which(df$Y == 1)
		roi.not.in = which(df$Y != 1)

		keepmat = expand.grid(low = c(0, 10, 20:30, 40, 45), 
			high = c(60, 70, 80:85, 100))
		keepmat$nvox_notroi = NA
		keepmat$nvox_roi = NA


		for (istop in seq(nrow(keepmat))){
			low = keepmat$low[istop]
			high = keepmat$high[istop]
			keepmat$nvox_roi[istop] = 
				sum( df$value[roi.in] >= low & 
					df$value[roi.in] <= high)
			keepmat$nvox_notroi[istop] = 
				sum( df$value[roi.not.in] >= low & 
					df$value[roi.not.in] <= high)				
			# cat(paste0(istop, " out of ", nrow(keepmat), "\n"))
		}
		nvox_notroi = length(roi.not.in)
		nvox_roi = length(roi.in)
		nvox = nvox_roi + nvox_notroi

		keepmat$pct_notroi = keepmat$nvox_notroi / nvox_notroi
		keepmat$pct_roi = keepmat$nvox_roi / nvox_roi
		# g = ggplot(keepmat, aes(x = pct_roi, y= pct_notroi, 
		# 	colour = factor(low), 
		# 	size = high)) + geom_point() 

		which.file = "ROI_Histograms.R"
		save(keepmat, nvox_notroi, nvox_roi, which.file, nvox,
			file = file.path(id.outdir, 
			paste0(img.stub, "_roi_keep_matrix", adder, ".Rda"))
			)

		######################################
		# Keep all ROI = 1, even if not inmask
		######################################	
		df = df[roi.in,]

		hval = hist(df$value, main=img.stub, breaks=200)
		dval = density(df$value)
		outname = nii.stub(basename(fdf$img[get.pred]))
		outname = file.path(id.outdir, 
			paste0(outname, "_ROI_values", adder, ".Rda"))
		vals = df$value
		save(vals, dval, hval, file=outname)	
		print(get.pred)
	}
	dev.off()

# }