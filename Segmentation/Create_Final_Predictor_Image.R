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
library(extrantsr)
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
correct = "Rigid"
# options = c("none", "N3", "N4", "N3_SS", "N4_SS",
#         "SyN", "SyN_sinc", "Rigid", "Affine", "Rigid_sinc", 
#         "Affine_sinc")
# options = c("none", "N3_SS", "N4_SS",
# 	"Rigid",  "Rigid_sinc")
options = "Rigid"

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

types = c("_zval2")
# , "_zval2"
# "_include_all", 
type = types[1]


keep.obj = ls()


	
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

	outfile = file.path(outdir, "111_Filenames_with_volumes.Rda")
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
	mod.filename = file.path(outdir, 
		paste0("Collapsed_Models", adder, ".Rda"))
	load(mod.filename)
	cns = colnames(res)

	get.pred <- as.numeric(Sys.getenv("SGE_TASK_ID"))
	if (is.na(get.pred)) get.pred = 16
	runpreds = seq(nrow(fdf))

	mycols = c("dice", "jaccard", "sens", "spec", "accur", 
		"truevol", "estvol")
	sim.res = matrix(NA, nrow=nrow(fdf), ncol = length(mycols))
	colnames(sim.res) = mycols

for (get.pred in runpreds){

	x = fdf[get.pred,]
	roi.fname = switch(correct,
		"none" = x$roi,
		"N3" = x$roi,
		"N4" = x$roi,
		"N3_SS" = x$roi,
		"N4_SS" = x$roi,
		"SyN" = x$synssroi,
		"SyN_sinc" = x$sinc_synssroi,
		"Rigid" = x$rig_ssroi,
		"Affine" = x$aff_ssroi,
		"Rigid_sinc" = x$sinc_rig_ssroi,
		"Affine_sinc" = x$sinc_aff_ssroi
		)		
	# print(get.pred)
	fname = switch(correct,
		"none" = x$img,
		"N3" = x$n3img,
		"N4" = x$n4img,
		"N3_SS" = x$n3ssimg,
		"N4_SS" = x$n4ssimg,
		"SyN" = x$synssimg,
		"SyN_sinc" = x$sinc_synssimg,
		"Rigid" = x$rig_ssimg,
		"Affine" = x$aff_ssimg,
		"Rigid_sinc" = x$sinc_rig_ssimg,
		"Affine_sinc" = x$sinc_aff_ssimg		
		)


	iddir = fdf$iddir[get.pred]
	id.outdir = fdf$outdir[get.pred]
	
	cut.filename = file.path(outdir, 
	paste0("Model_Cutoffs", adder, ".Rda"))

	load(file=cut.filename)

	scut.filename = file.path(outdir, 
	paste0("Smooth_Model_Cutoffs", adder, ".Rda"))

	load(file=scut.filename)

	names(all.scuts) = cns

	# id.outdir = x$outdir
	# predname = nii.stub(basename(x$img))
	# predname = file.path(id.outdir, 
	# 	paste0(predname, "_predictors", adder, ".Rda"))
	# load(predname)
	# df = img.pred$df
	# rm(list="img.pred")

	# for (ipred in seq(ncol(preds))){
	cn = "mod_agg"

		roi = readNIfTI(roi.fname, reorient=FALSE)
		img = readNIfTI(fname, reorient=FALSE)
		smimg = nii.stub(fdf$img[get.pred], bn=TRUE)
		smimg = file.path(id.outdir, 
			paste0(smimg, "_", cn, adder))
		smimg = paste0(smimg, "_smoothed", adder)
		sm.img = readNIfTI(fname = smimg, reorient=FALSE)
		sm.img = niftiarr(sm.img, sm.img > all.scuts[cn])
		outimg = nii.stub(fdf$img[get.pred], bn=TRUE)
		outimg = file.path(id.outdir, 
			paste0(outimg, "_", cn, adder, "_prediction"))
		# outimg = paste0(outimg, ".nii.gz")
		if (file.exists(paste0(outimg, ".nii.gz.nii.gz"))){
			file.remove(paste0(outimg, ".nii.gz.nii.gz"))
		}
		writeNIfTI(sm.img, filename = outimg )

		myres = extrantsr:::sim(dman = c(roi), 
			dauto = c(sm.img), 
			dim.dman = dim(roi),
			dim.dauto = dim(sm.img))
		myres = unlist(myres)
		sim.res[get.pred, mycols] = myres[mycols]

	# }
	print(get.pred)
}

sim.res = cbind(fdf[runpreds,], sim.res)

outfile = file.path(outdir, 
	paste0("Overlap_Measures", adder, "_", cn, ".Rda"))

save(sim.res, file=outfile)

# }
# }

# res = t(res)
# sres = t(sres)



