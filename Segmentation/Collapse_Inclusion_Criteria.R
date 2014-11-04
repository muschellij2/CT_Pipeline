###################################################################
## This code is for inclusion criteria colapsing
##
## Author: John Muschelli
## Last updated: May 20, 2014
###################################################################
###################################################################
rm(list=ls())
library(plyr)
library(ggplot2)
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
# options = c("none", "N3", "N4", "N3_SS", "N4_SS",
#         "SyN", "SyN_sinc", "Rigid", "Affine", "Rigid_sinc", 
#         "Affine_sinc")
options = c("none", "N3_SS", "N4_SS", 
		"Rigid", "Rigid_sinc")

icorr <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(icorr)) icorr = 1
correct = options[icorr]

keep.obj = ls()

l = list()
for (icorr in seq_along(options)){
	correct = options[icorr]
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


	predname = file.path(outdir, 
		paste0("Aggregate_cutoff_results", adder, ".Rda"))

	load(file = predname)

	cmeans = colMeans(ffdf[, cn])
	cmaxs = colMaxs(ffdf[, cn])
	cmins = colMins(ffdf[, cn])

	x = data.frame(mean=t(t(cmeans)), max=t(t(cmaxs)), 
		min = t(t(cmins)))
	x$var=  rownames(x)
	x$cut = gsub("(.*)[.](out|reduced)", "\\1", x$var)
	x$type = gsub("(.*)[.](out|reduced)", "\\2", x$var)
	x$cut = gsub(".cutoff", "", x$cut)
	rownames(x) = NULL
	x$var = NULL
	res = reshape(x, direction = "wide", idvar = c("cut"), 
		timevar = "type")

	l[[icorr]] = list(res=res, ffdf = ffdf)
	# l[[icorr]] = list(ffdf = ffdf)
}

x = l[[1]]$res
plot( max.out ~ min.reduced, data=x)
text(x= x$min.reduced, y = x$max.out, labels = x$cut)

plot( mean.out ~ mean.reduced, data=x)
text(x= x$mean.reduced, y = x$mean.out, labels = x$cut)