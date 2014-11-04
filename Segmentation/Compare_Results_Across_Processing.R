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
        "SyN", "SyN_sinc", "Rigid", "Affine", "Rigid_sinc", 
        "Affine_sinc")
nopts = length(options)


types = c("", "_include", "_zval", "_zval2")
type = types[1]

#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )

outfile = file.path(outdir, "111_Filenames.Rda")
load(file = outfile)

mod.filename = file.path(outdir, 
	paste0("Collapsed_Models.Rda"))
x = load(mod.filename)	
cn = colnames(res)

nopred = seq(non.aggmods)
vd = vol.data
vd = vd[-nopred,]
nr = nrow(vd)
valid.ind = ceiling(nr/2)
test.ind = seq( valid.ind +1, nr)
valid.ind = seq(1, valid.ind)
group = "Training"

ffdf = fdf[-nopred, ]
if (group == "Training"){
	subset.ind = valid.ind
}
if (group == "Test"){
	subset.ind = test.ind
}

subset.ids = ffdf$id[subset.ind]

rm(list=c(x, "vd"))

p = length(cn)

m.pauc = matrix(NA, nrow=length(options), ncol = p)
colnames(m.pauc) = cn

m.spauc = m.pauc
m.vdiff = m.pauc
m.svdiff = m.vdiff

truevols = matrix(NA, ncol=length(options), 
	nrow = length(subset.ind))
colnames(truevols) = options
icorr = 1




for (icorr in seq(nopts)){
	
	correct = options[icorr]
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
		"Affine" = "_Affine",
		"Rigid_sinc" = "_Rigid_sinc",
		"Affine_sinc" = "_Affine_sinc")



	mod.filename = file.path(outdir, 
		paste0("Collapsed_Models", adder, ".Rda"))
	load(mod.filename)	

	outfile = file.path(outdir, 
		paste0("Model_performance_results", adder, type, 
			".Rda")
		)
	x=load(file = outfile)

	vsd = vol.sdata
	vsd = vsd[-nopred,]
	vsd = vsd[, !(colnames(vsd) %in% "truth")]

	###############################
	# Separate data into validation and test
	###############################
	vd = vol.data
	vd = vd[-nopred,]
	nr = nrow(vd)	
	truevol = vd[,"truth"]
	vd = vd[, !(colnames(vd) %in% "truth")]
	voldiff = vd - truevol
	adiff = abs(voldiff)

	######################################
	# Taking volume differences from truth and subsetting
	######################################
	cut.diff = cut.vol.datas[-nopred, 
			!colnames(cut.vol.datas) %in% "truth"] - 
		truevol
	cut.adiff = abs(cut.diff)

	cut.sdiff = cut.vol.sdatas[-nopred, 
			!colnames(cut.vol.sdatas) %in% "truth"] -
		truevol
	cut.sadiff = abs(cut.sdiff)

	cut.tsdiff = cut.vol.tsdatas[-nopred, 
			!colnames(cut.vol.sdatas) %in% "truth"] -
		truevol
	cut.tsadiff = abs(cut.tsdiff)				

	svoldiff = vsd - truevol
	sadiff = abs(svoldiff) 

	######################################
	# The ones with a's are absolute, should just do in plot code
	######################################
	vol = voldiff[subset.ind,]

	avol = adiff[subset.ind,]

	svol = svoldiff[subset.ind,]
	savol = sadiff[subset.ind,]

	cut.svol = cut.sdiff[subset.ind,]
	cut.savol = cut.sadiff[subset.ind,]

	cut.valid.tsvol = cut.tsdiff[subset.ind,]
	cut.valid.tsavol = cut.tsadiff[subset.ind,]

	cut.valid.vol = cut.diff[subset.ind,]
	cut.valid.avol = cut.adiff[subset.ind,]

	######################################
	# Getting pAUC
	######################################
	res = res[-nopred, ]
	sres = sres[-nopred, ]

	res = res[subset.ind,]
	sres = sres[subset.ind,]

	######################################
	# Getting Accuracy
	######################################
	all.accs = all.accs[-nopred, ]
	all.saccs = all.saccs[-nopred, ]

	acc = all.accs[subset.ind,]
	sacc = all.saccs[subset.ind,]

	cut.all.accs = mod.accs[-nopred, ]
	cut.all.saccs = mod.saccs[-nopred, ]

	cut.acc = cut.all.accs[subset.ind,]
	cut.sacc = cut.all.saccs[subset.ind,]


	######################################
	# Getting Benchmarks
	######################################
	benches = benches[-nopred]
	bench = benches[subset.ind]

	bench.m = matrix(bench, nrow=length(bench), 
		ncol = ncol(acc), byrow=FALSE)


	fcn = colMedians
	m.pauc[icorr, ] = do.call(fcn, list(res))
	m.spauc[icorr, ] = do.call(fcn, list(sres))

	m.vdiff[icorr, ] = do.call(fcn, list(avol))
	m.svdiff[icorr, ] = do.call(fcn, list(savol))
	truevols[,icorr] = truevol[subset.ind]
	print(icorr)
}

tvol = truevols[, "none"]
tvols = truevols[, colnames(truevols) %in% c("none")]
vdiff = (tvols - tvol)
vdat = melt(vdiff)

which(m.pauc == max(m.pauc), arr.ind=TRUE)
max(m.pauc)
which(m.spauc == max(m.spauc), arr.ind=TRUE)
max(m.spauc)
which(m.vdiff == min(m.vdiff), arr.ind=TRUE)
min(m.vdiff)
which(m.svdiff == min(m.svdiff), arr.ind=TRUE)
min(m.svdiff)
# }