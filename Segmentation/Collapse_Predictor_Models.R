###########################################################
## This code is for collapsing predictor models
##
## Author: John Muschelli
## Last updated: May 20, 2014
###########################################################
###########################################################
rm(list=ls())
library(plyr)
library(cttools)
library(fslr)
library(ROCR)
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

correct = "Rigid"
# options = c("none", "N3", "N4", "N3_SS", "N4_SS",
#         "SyN", "SyN_sinc", "Rigid", "Affine", 
# "Rigid_sinc", 
#         "Affine_sinc")
options = c("none", "N3_SS", "N4_SS", 
	"Rigid", "Rigid_sinc")

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
		"Affine" = "_Affine",
		"Rigid_sinc" = "_Rigid_sinc",
		"Affine_sinc" = "_Affine_sinc")

    filename = file.path(outdir, 
        paste0("Result_Formats", adder, ".Rda"))
    xxx = load(filename)


	all.dat = llply(seq(lmod), function(imod){
			mod.outdir = fdf.run$outdir[imod]
			moddname = nii.stub(basename(fdf.run$img[imod]))
			moddname = file.path(mod.outdir, 
				paste0(moddname, "_models", adder, ".Rda"))
			print(moddname)
			x = load(moddname)
			mod = mods$mod
			if (is.vector(acc)) {
				acc = t(as.matrix(acc))
			}
			if (is.vector(pauc.cut)) {
				pauc.cut = t(as.matrix(pauc.cut))
			}			
			cutoff = acc[, 'cutoff']
			pauc.cutoff = pauc.cut[, "cutoff"]
			sens.cutoff = sens.cut[, "cutoff"]
			dice.cutoff = dice.coef[, "cutoff"]
			# gam.cutoff = gam.acc[, "cutoff"]
			# gam.pauc.cutoff = gam.pauc.cut[, "cutoff"]

			l = list(mod= mod, cutoff= cutoff, 
				pauc.cutoff = pauc.cutoff,
				sens.cutoff = sens.cutoff,
				dice.cutoff = dice.cutoff
				# gam.cutoff = gam.cutoff,
				# gam.pauc.cutoff = gam.pauc.cutoff,
				# gam.mod = gam.mod
				)
			return(l)
	}, .progress = "text")

	all.mods = llply(all.dat, function(x){
		x$mod
	})

	imod = lmod
	mod.outdir = fdf.run$outdir[imod]
	moddname = nii.stub(basename(fdf.run$img[imod]))
	moddname = file.path(mod.outdir, 
		paste0(moddname, "_models", adder, ".Rda"))
	print(moddname)	
	xx = load(moddname)


	# all.gam.mods = llply(all.dat, function(x){
	# 	x$gam.mod
	# })


	all.cutoffs = laply(all.dat, `[[`, "cutoff")
	all.pauc.cutoffs = laply(all.dat, `[[`, "pauc.cutoff")
	all.sens.cutoffs = laply(all.dat, `[[`, "sens.cutoff")
	all.dice.cutoffs = laply(all.dat, `[[`, "dice.cutoff")

	names(all.mods) = names(all.cutoffs) = 
	names(all.pauc.cutoffs) = names(all.sens.cutoffs) = 
		colnames(res)[seq(lmod)]
	names(all.dice.cutoffs) = names(all.sens.cutoffs) 
	#########
	# need to get cutoffs too
	##########


	filename = file.path(outdir, 
		paste0("Collapsed_Models", adder, ".Rda"))
	save(all.mods, res, vol.data, vol.sdata, sres, fdf.run,
		runpreds, run.ind,
		all.cutoffs, 
		all.pauc.cutoffs,
		all.dice.cutoffs, 
		all.sens.cutoffs,
		gam.mod, 
		gam.acc,  
		gam.pauc, 
        gam.sens.cut,
        gam.dice.coef,	
        gam.pauc.cut,
        rf.sens.cut, 
        rf.dice.coef,
        rf.mod, rf.acc, rf.pauc, 
        rf.pauc.cut,        
		lmod, non.aggmods, file=filename)
	print(correct)
# }