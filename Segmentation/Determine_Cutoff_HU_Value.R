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
        "SyN", "SyN_sinc", "Rigid", "Affine", "Rigid_sinc", 
        "Affine_sinc")
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


	mod.filename = file.path(outdir, 
		paste0("Collapsed_Models", adder, ".Rda"))
	load(mod.filename)
	nopred = seq(non.aggmods)
	
	fdf.run = fdf[-nopred,]
	nr = nrow(fdf)
	valid.ind = ceiling(nr/2)
	test.ind = seq( valid.ind +1, nr)
	valid.ind = seq(1, valid.ind)

	mygroup = "Validation"
	fdf = fdf.run[valid.ind, ]


	##############################
	# Keeping files where predictors exist
	##############################
	outfiles = nii.stub(basename(fdf$img))
	outfiles = paste0(outfiles, "_roi_keep_matrix", adder, ".Rda")
	outfiles = file.path(fdf$outdir, outfiles)
	stopifnot(file.exists(outfiles))

	get.pred = 1
	x = fdf[get.pred,]

	all.mat = NULL

	for (get.pred in seq(nrow(fdf))){

		iddir = fdf$iddir[get.pred]
		id.outdir = fdf$outdir[get.pred]
		img.stub = nii.stub(fdf$img[get.pred], bn=TRUE)
		rdaname = file.path(id.outdir, 
			paste0(img.stub, "_roi_keep_matrix", adder, ".Rda"))
		load(rdaname)
		keepmat$nvox = nvox
		keepmat$img = img.stub
		keepmat$total_roi = nvox_roi
		keepmat$total_notroi = nvox_notroi
		all.mat = rbind(all.mat, keepmat)
		print(get.pred)
	}

	agg = ddply(all.mat, .(low, high), summarise,
		nvox_notroi = sum(nvox_notroi), 
		nvox_roi = sum(nvox_roi), 
		total_notroi = sum(total_notroi), 
		total_roi = sum(total_roi), 		
		nvox = sum(nvox)
		)



	agg$pct_roi = agg$nvox_roi / agg$total_roi
	agg$pct_notroi = agg$nvox_notroi / agg$total_notroi

	agg.good = agg[ agg$high > 70, ]

	aa = agg.good[agg.good$low == 30,]

	pdfname = file.path(outdir, 
	    paste0("HU_Threshold", adder, "_", mygroup, ".pdf"))	
	pdf(pdfname)
		g = ggplot(agg, aes(x = pct_roi, y= pct_notroi, 
		colour = factor(low), 
		size = high)) + geom_point() 
		g.good = g %+% agg.good
		ga = g %+% aa
		print(g)
		print(g.good)
		print(ga)
	dev.off()
}

	# g = ggplot(all.mat, aes(x = pct_roi, y= pct_notroi, 
	# 	colour = factor(low), 
	# 	size = high)) + geom_point() 


	# g = ggplot(all.mat, aes(x = pct_roi, y= pct_notroi)) + 
	# 	facet_wrap(low ~ high) + geom_point() 	

	# pdfname = file.path(outdir, 
	#     paste0("ROI_Histogram", adder, ".pdf"))
	# pdf(pdfname)	
	# dev.off()

# }