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

types = c("_zval2", "_zval_all", "_zval2_medztemp")
# , "_zval2"ls
# "_include_all", 
type = types[1]


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

	outfile = file.path(outdir, 
		"111_Filenames_with_volumes_stats.Rda")
	load(file = outfile)


	fdf$imgfname = switch(correct, 
		"none"= fdf$ssimg,
		"N3"=fdf$n3img,
		"N4" = fdf$n4img,
		"N3_SS" = fdf$n3ssimg,
		"N4_SS" = fdf$n4ssimg, 
		"SyN" = fdf$synssimg,
		"SyN_sinc" = fdf$sinc_synssimg,
		"Rigid" = fdf$rig_ssimg,
		"Affine" = fdf$aff_ssimg,
		"Rigid_sinc" = fdf$sinc_rig_ssimg,
		"Affine_sinc" = fdf$sinc_aff_ssimg)		

	fdf$roifname = switch(correct, 
		"none"= fdf$roi,
		"N3"=fdf$roi,
		"N4" = fdf$roi,
		"N3_SS" = fdf$roi,
		"N4_SS" = fdf$roi, 
		"SyN" = fdf$synssroi,
		"SyN_sinc" = fdf$sinc_synssroi,
		"Rigid" = fdf$rig_ssroi,
		"Affine" = fdf$aff_ssroi,
		"Rigid_sinc" = fdf$sinc_rig_ssroi,
		"Affine_sinc" = fdf$sinc_aff_ssroi)		

	fdf$preddir = switch(correct, 
		"none"= fdf$roi,
		"N3"=fdf$roi,
		"N4" = fdf$roi,
		"N3_SS" = fdf$roi,
		"N4_SS" = fdf$roi, 
		"SyN" = fdf$synssroi,
		"SyN_sinc" = fdf$sinc_synssroi,
		"Rigid" = fdf$rig_ssroi,
		"Affine" = fdf$aff_ssroi,
		"Rigid_sinc" = fdf$sinc_rig_ssroi,
		"Affine_sinc" = fdf$sinc_aff_ssroi)	



	makedir = sapply( fdf$outdir, function(x) {
		if (!file.exists(x)){
			dir.create(x, showWarnings =FALSE)
		}
	})
	irow = 1
	x = fdf[irow,]


	# load(file = file.path(outdir, "Segmentation_Models.Rda"))
	##############################
	# Run lmod number of models - not all the models - leave out
	##############################
#"223-369"
	get.pred <- as.numeric(Sys.getenv("SGE_TASK_ID"))
	if (is.na(get.pred)) get.pred = 71
	x = fdf[get.pred,]
	print(get.pred)

# for (get.pred in runpreds){

	iddir = x$iddir
	id.outdir = x$outdir
	predname = nii.stub(basename(x$img))
	predname = file.path(id.outdir, 
		paste0(predname, "_predictors", adder, ".Rda"))
	load(predname)
	df = img.pred$df
	stopifnot(all(df$Y %in% c(0, 1)))

	cut.filename = file.path(outdir, 
	paste0("Model_Cutoffs", adder, ".Rda"))

	xcut = load(file=cut.filename)

	scut.filename = file.path(outdir, 
	paste0("Smooth_Model_Cutoffs", adder, ".Rda"))

	xscut = load(file=scut.filename)

	cnames = names(all.cuts)

	img = readNIfTI(x$imgfname, reorient = FALSE)
	roi = readNIfTI(x$roifname, reorient = FALSE)
	roi = cal_img (roi > 0)
	xyz = cog(roi, ceil=TRUE)


	# for (type in types){
	# typeimg = file.path(x$outdir, 
 #    	paste0(nii.stub(x$imgfname, bn=TRUE), 
 #    		"_zval2_medztemp"))
	# timg = readNIfTI(typeimg, reorient=FALSE)	    

		ipred = which(cnames == "mod_agg")
		################################
		# Running with variable slice thickness
		################################
		# for (ipred in seq_along(cnames)){


		if (correct %in% c("none", "N3_SS", "N4_SS")) {
			ofile = file.path( x$outdir, 
			paste0(nii.stub(x$ss, bn=TRUE), 
				"_template_zscore.nii.gz" ))
		}
		if (correct %in% c("Rigid", "Rigid_sinc")) {
		    outputdir = paste0("Rigid_Registered")
		    outputdir = file.path(x$iddir, outputdir)
		    ofile = file.path( outputdir, 
		        paste0(nii.stub(x$ss, bn=TRUE), 
		            "_template_zscore_", correct ))
		}

			cn = cnames[ipred]
			# pred.imgs[[ipred]] = img
			xout = outimg = nii.stub(x$img, bn=TRUE)
			outimg = file.path(id.outdir, 
				paste0(outimg, "_", cn, adder))	

			########################
			# Read in Z-score template image
			########################
			zimg = readNIfTI(ofile, reorient=FALSE)
			rzimg = robust_window(zimg)

			reg.img = readNIfTI(outimg, reorient=FALSE)
			reg.thresh = all.cuts[cn]
			reg.mask = cal_img(reg.img > reg.thresh)

			total.mask = cal_img(reg.img > 0)


			outimg = paste0(outimg, "_smoothed", adder, type)
			sm.img = readNIfTI(outimg, reorient=FALSE)

			thresh = all.scuts[cn]
			sm.mask = cal_img(sm.img > thresh)

			mask.overlay(img, sm.mask, window=c(0, 100),
				xyz=xyz)
			# x11()
			mask.overlay(img, reg.mask, window=c(0, 100),
				xyz=xyz)
			# x11()
			mask.overlay(img, roi, window=c(0, 100),
				xyz=xyz)

			diffimg = niftiarr(roi, (sm.mask - roi))
			diffimg[diffimg == -1] = 2
			mask.overlay(img, diffimg, window=c(0, 100), 
				col.y = c("blue", "red"), 
				xyz=cog(diffimg, ceil=TRUE), 
				ybreaks=c(0, 1, 2),
				zlim.y = c(0,2))
		

			fname= nii.stub(basename(x$img))
			fname = paste0(fname, "_predictors", 
				adder, "_training.Rda")
			outfile = file.path(x$outdir[get.pred], fname)
			onames = load(file=outfile)			
		# }
		# close(pb)




