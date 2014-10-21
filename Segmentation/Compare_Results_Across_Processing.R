##################################################################
## This code is comparing results (AUC and such) across procs
##
## Author: John Muschelli
## Last updated: May 20, 2014
##################################################################
##################################################################
rm(list=ls())
library(plyr)
library(cttools)
library(ROCR)
library(matrixStats)
library(reshape2)
library(ggplot2)
library(fslr)
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

outdir = file.path(basedir, "results")
correct = "none"
options = c("none", "N3", "N4", "N3_SS", "N4_SS",
		"SyN", "SyN_sinc", "Rigid", "Affine")

# for (correct in options){
	
	print(correct)
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


	mod.filename = file.path(outdir, 
		paste0("Collapsed_Models", adder, ".Rda"))
	load(mod.filename)	

	outfile = file.path(outdir, 
		paste0("Model_performance_results", adder, ".Rda")
		)
	load(file = outfile)


	nopred = seq(non.aggmods)
	vd = vol.data
	vd = vd[-nopred,]
	nr = nrow(vd)
	valid.ind = ceiling(nr/2)
	test.ind = seq( valid.ind +1, nr)
	valid.ind = seq(1, valid.ind)


	ffdf = fdf[-nopred, ]
	group = "Training"
	if (group == "Training"){
		subset.ind = valid.ind
	}
	if (group == "Test"){
		subset.ind = test.ind
	}

	vsd = vol.sdata
	vsd = vsd[-nopred,]
	vsd = vsd[, !(colnames(vsd) %in% "truth")]

	###############################
	# Separate data into validation and test
	###############################
	truevol = vd[,"truth"]
	vd = vd[, !(colnames(vd) %in% "truth")]
	voldiff = vd - truevol
	adiff = abs(voldiff)

	cut.adiff = abs(
		cut.vol.datas[-nopred, 
			!colnames(cut.vol.datas) %in% "truth"] - 
		truevol
		)
	cut.sadiff = abs(
		cut.vol.sdatas[-nopred, 
			!colnames(cut.vol.sdatas) %in% "truth"] -
		truevol
		)	

	sadiff = abs(vsd - truevol) 
	svol = sadiff[subset.ind,]

	vol = adiff[subset.ind,]

	cut.svol = cut.sadiff[subset.ind,]

	cut.valid.vol = cut.adiff[subset.ind,]

	res = res[-nopred, ]
	sres = sres[-nopred, ]

	all.accs = all.accs[-nopred, ]
	all.saccs = all.saccs[-nopred, ]
	
	benches = benches[-nopred]

	res = res[subset.ind,]

	sres = sres[subset.ind,]

	acc = all.accs[subset.ind,]
	sacc = all.saccs[subset.ind,]

	cut.all.accs = mod.accs[-nopred, ]
	cut.all.saccs = mod.saccs[-nopred, ]	

	cut.acc = cut.all.accs[subset.ind,]
	cut.sacc = cut.all.saccs[subset.ind,]


	bench = benches[subset.ind]

	bench.diff = acc - bench
	n_above = colSums(bench.diff > 0)

	above_bench = apply(bench.diff > 0, 2, all)
	best.acc = which(above_bench)
	biggest.acc_diff = which.max(colMeans(bench.diff))

	best.res= apply(res, 1, which.max)
	best.sres= apply(sres, 1, which.max)

	# sacc = valid.sacc
	# acc = valid.acc
	# inds = valid.ind
	# benches = valid.bench
	# vol = valid.vol
	# svol = valid.svol
	# sres = valid.sres

	subset.ids = ffdf$id[subset.ind]



	

# }