
##########################################################
## This code is for final prediction volumes of 
# Image Segmentation of CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
##########################################################
##########################################################
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



###########################
# Varying all the parameters
###########################
# types = c("_zval2", '_zval2_medztemp')
types = "_zval2"
# "_include_all", 
type = types[1]
		
###########################
# Varying the cutoff values
###########################
cutter = c("", "_dice")
icut = cutter[1]

###########################
# Varying the models
###########################
cns = c("mod_agg", "gam", "rf")
# cns = "gam"
cn = cns[1]

###########################
# making scenarios
###########################
options = c("Rigid", "none")

scenarios = expand.grid(icut = cutter, 
	type = types,
	cn = cns, 
	correct = options, 
	stringsAsFactors = FALSE)

scenarios = scenarios[ scenarios$type == "_zval2", ]

# 

iscen <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(iscen)) iscen = 1

cn = scenarios$cn[iscen]
type = scenarios$type[iscen]
icut = scenarios$icut[iscen]
correct = scenarios$correct[iscen]

# scenarios$insuf = paste0("_", scenarios$cn, adder, 
# 						"_smoothed", adder, scenarios$type)
# scenarios$outsuf = paste0("_", scenarios$cn, adder, 
# 	scenarios$type, 
# 	scenarios$icut, "_prediction")




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




keep.obj = ls()


mod.filename = file.path(outdir, 
	paste0("Collapsed_Models", adder, ".Rda"))
load(mod.filename)
cns = colnames(res)

	
cut.filename = file.path(outdir, 
paste0("Model_Cutoffs", adder, ".Rda"))

load(file=cut.filename)

scut.filename = file.path(outdir, 
paste0("Smooth_Model_Cutoffs", adder, ".Rda"))

load(file=scut.filename)
names(all.sdice.cuts) = names(all.scuts) = cns

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

# for (iscen in seq(nrow(scenarios))){


	get.pred = 38
	runpreds = seq(nrow(fdf))

	mycols = c("dice", "jaccard", "sens", "spec", "accur", 
		"truevol", "estvol")

	sim.res = matrix(NA, nrow=nrow(fdf), 
		ncol = length(mycols))
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


		if (icut == ""){
			cutvec = all.scuts
		}
		if (icut == "_dice"){
			cutvec = all.sdice.cuts
		}

		# for (ipred in seq(ncol(preds))){

		####################################
		# Need to read in training results for best.mat 
		# and 
		# get a predictor for each to calculate volumes
		####################################

		roi = readNIfTI(roi.fname, reorient=FALSE)
		roi = cal_img(roi > 0.5)


		# img = readNIfTI(fname, reorient=FALSE)
		smimg = nii.stub(fdf$img[get.pred], 
			bn=TRUE)
		smimg = file.path(id.outdir, 
			paste0(smimg, "_", cn, adder, 
				"_smoothed", adder, type))
		sm.img = readNIfTI(fname = smimg, 
			reorient=FALSE)
		sm.img = niftiarr(sm.img, 
			sm.img > cutvec[cn])

		outimg = nii.stub(fdf$img[get.pred], bn=TRUE)
		outimg = file.path(id.outdir, 
			paste0(outimg, "_", cn, adder, type, icut, 
				"_prediction"))

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

		print(get.pred)
	}
	# }
	sim.res = cbind(fdf[runpreds,], sim.res)

	outfile = file.path(outdir, 
		paste0("Overlap_Measures", adder, "_", cn, 
			type, icut,
			".Rda"))

	save(sim.res, file=outfile)

# }
# }

# res = t(res)
# sres = t(sres)



