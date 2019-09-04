###################################################################
## This code is for unsmoothed cutoffs for models
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
library(mgcv)
library(extrantsr)
library(randomForest)
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

segdir = file.path(progdir, "Segmentation")
source(file.path(segdir, "performance_functions.R"))

correct = "Rigid_sinc"
# options = c("none", "N3", "N4", "N3_SS", "N4_SS",
#         "SyN", "SyN_sinc", "Rigid", "Affine", "Rigid_sinc", 
#         "Affine_sinc")
options = c("none",  "N3_SS", "N4_SS",
    "Rigid", "Rigid_sinc")


#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )

outfile = file.path(outdir, "111_Filenames_with_volumes_stats.Rda")
load(file = outfile)

icorr <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(icorr)) icorr = 1
correct = options[icorr]

# for (correct in options){
    if ("all.df" %in% ls()){
        rm(list="all.df")
    }
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


    ##############################
    # Keeping files where predictors exist
    ##############################
    outfiles = nii.stub(basename(fdf$img))
    outfiles = paste0(outfiles, "_predictors", adder, ".Rda")
    outfiles = file.path(fdf$outdir, outfiles)
    fdf = fdf[file.exists(outfiles), ]

    # load(file = file.path(outdir, "Segmentation_Models.Rda"))
    ##############################
    # Run lmod number of models - not all the models - leave out
    ##############################
	mod.filename = file.path(outdir, 
		paste0("Collapsed_Models", adder, ".Rda"))
	x= load(mod.filename)


    fdf.run = fdf[run.ind, ]

    moddname = nii.stub(basename(fdf.run$img))
    moddname = file.path(fdf.run$outdir, 
        paste0(moddname, "_predictors", adder, ".Rda"))


    fname = file.path(outdir, 
        paste0("Candidate_Aggregate_data", adder, ".Rda"))
    load(fname)

    test$mode = fdf$mode[ match(test.img, fdf$img)]

	preds = t(laply(all.mods, short_predict, 
        newdata= test,
        .progress = "text"))

    rowMean = rowMeans(preds)
    rowMed = rowMedians(preds)
    rowMax = rowMaxs(preds)
    rowMin = rowMins(preds)
    rowProd = exp(rowSums(log(preds)))
    nr = nrow(preds)
    rowGeom = rowProd ^ (1/nr) 
    # rowGeom = exp(rowMeans(log(preds)))
    rf = predict(rf.mod, 
        newdata = test, 
        type="prob")[, "1"]

    preds = cbind(rowMean, rowMed, 
        rowMax, rowMin,
        rowProd, rowGeom, rf)
    rownames(preds) = NULL


	#### only values 0 to 100 are likely to be ROI
	# preds = preds * df$in0100
	preds = preds * test.mult.df$multiplier
	fpr.stop = 0.01

	r.res = alply(preds, 2, function(nd){
		pred <- prediction( nd, test$Y)

        N = nrow(test)
        sens.cut = get_senscut(pred, fpr.stop=fpr.stop) 
            # N = N,
            # predictions=nd, Y = test$Y)
        pauc = get_pauc(pred, fpr.stop=fpr.stop)
        dice.coef = get_max_dice(pred)
        acc = get_acc(pred)

		perf <- performance(pred, "tpr", "fpr")
		pauc.cut = t(opt.cut(perf, pred))

		list(acc = acc, pauc = pauc, pauc.cut = pauc.cut,
            sens.cut = sens.cut, 
            dice.coef = dice.coef)
	}, .progress = "text")

	accs = sapply(r.res, `[[`, "acc")
	paucs = sapply(r.res, `[[`, "pauc")
	pauc.cuts = sapply(r.res, `[[`, "pauc.cut")[3,]
    sens.cuts = sapply(r.res, `[[`, "sens.cut")[4,]
    dice.cuts = sapply(r.res, `[[`, "dice.coef")[2,]

	colnames(accs) = names(paucs) = names(pauc.cuts) = 
    names(dice.cuts) = names(sens.cuts) = colnames(preds)

	cuts = accs[1, ]
	all.cuts = c(all.cutoffs, cuts, gam=gam.acc[, "cutoff"])
	all.pauc.cuts = c(all.pauc.cutoffs, pauc.cuts, 
        gam=gam.pauc.cut[, "cutoff"])
    all.pauc.cuts = c(all.pauc.cutoffs, pauc.cuts, 
        gam=gam.pauc.cut[, "cutoff"])   
    all.sens.cuts = c(all.sens.cutoffs, sens.cuts, 
        gam=gam.sens.cut[, "cutoff"])            
    all.dice.cuts = c(all.dice.cutoffs, dice.cuts, 
        gam=gam.dice.coef[, "cutoff"])                
    names(all.pauc.cuts)[length(all.pauc.cuts)] = "gam"
    names(all.sens.cuts)[length(all.sens.cuts)] = "gam"
    names(all.dice.cuts)[length(all.dice.cuts)] = "gam"
    names(all.cuts)[length(all.cuts)] = "gam"

    all.cuts[is.infinite(all.cuts)] = 1
    all.pauc.cuts[is.infinite(all.pauc.cuts)] = 1
    all.sens.cuts[is.infinite(all.sens.cuts)] = 1
    all.dice.cuts[is.infinite(all.dice.cuts)] = 1

	mod.filename = file.path(outdir, 
		paste0("Model_Cutoffs", adder, ".Rda"))

	save(all.cuts, all.pauc.cuts, 
        all.sens.cuts,
        all.dice.cuts, 
        file=mod.filename)

# }