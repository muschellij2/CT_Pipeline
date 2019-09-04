###################################################################
## This code is for kmeans clustering
##
## Author: John Muschelli
## Last updated: May 20, 2014
###################################################################
###################################################################
rm(list=ls())
library(plyr)
library(cttools)
library(fslr)
library(matrixStats)
library(getopt)
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

segdir = file.path(progdir, "Segmentation")
source(file.path(segdir, "performance_functions.R"))

outdir = file.path(basedir, "results")

correct = "Rigid"

options = c("none", "N3_SS", "N4_SS",
	"Rigid",  "Rigid_sinc")
# options = c("none", "Rigid_sinc")
# options = "Rigid"

spec = matrix(c(
	'correct', 'c', 1, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

correct = opt$correct
print(opt)

my.tab <- function(
  x, 
  y, 
  dnames=c("x", "y")) {
  x = as.numeric(x)
  y = as.numeric(y)
  stopifnot(all(unique(c(x,y)) %in% c(0, 1, NA)))
  tt = sum(x * y)
  t1=sum(x)
  t2=sum(y)
  tab = matrix(c(length(x)-t1-t2+tt,  t1-tt, t2-tt, tt), 2, 2)
  n = list(c("FALSE", "TRUE"), c("FALSE", "TRUE"))
  names(n) = dnames
  dimnames(tab) = n
  tab = as.table(tab)
  return(tab) 
}


keep.obj = ls()

# for (correct in options){
	
	all.obj = ls()
	rm.obj = all.obj[!(all.obj %in% c(keep.obj, "keep.obj"))]
	rm(list=rm.obj)
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
	outfiles = nii.stub(basename(fdf$img))
	outfiles = paste0(outfiles, "_predictors", adder, ".Rda")
	outfiles = file.path(fdf$outdir, outfiles)
	stopifnot(file.exists(outfiles))

	get.pred <- as.numeric(Sys.getenv("SGE_TASK_ID"))
	if (is.na(get.pred)) get.pred = 71
	x = fdf[get.pred,]
	print(get.pred)

# for (get.pred in runpreds){

	iddir = fdf$iddir[get.pred]
	id.outdir = fdf$outdir[get.pred]
	predname = nii.stub(basename(fdf$img[get.pred]))
	predname = file.path(id.outdir, 
		paste0(predname, "_predictors", adder, ".Rda"))
	load(predname)
	df = img.pred$df
	stopifnot(all(df$Y %in% c(0, 1)))

	remake_pred = function(nim, keep.ind, vec){
		nim = niftiarr(nim, array(0, dim=dim(nim)))
		nim[keep.ind] = vec
		nim = cal_img(nim)
		nim
	}
	nim = img.pred$nim
	keep.ind = img.pred$keep.ind
	# rm(list="img.pred")
	d = df[ keep.ind, ]
	Y = d$Y
	d$Y = d$mask = NULL
	d$dist_centroid = NULL
	img = remake_pred(nim, keep.ind, d$value)
	roi = remake_pred(nim, keep.ind, Y)
 	xyz=cog(roi, ceil=TRUE)

	km = kmeans(d, centers = 4);
	clustnum = which.max(km$centers[, "thresh"])
	clust_img = remake_pred(nim, keep.ind, km$cluster == clustnum)
	mask.overlay(img, clust_img, window=c(0, 100), xyz=xyz)

	table(km$cluster == clustnum, Y)
# }
