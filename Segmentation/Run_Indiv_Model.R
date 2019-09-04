###################################################################
## This code is for Image Segmentation of CT
## The code is R but calls out FSL
##
## Author: John Muschelli
## Last updated: May 20, 2014
###################################################################
###################################################################
rm(list=ls())
library(plyr)
library(cttools)
library(devtools)
library(ROCR)
library(fslr)
library(mgcv)
library(extrantsr)
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
outdir = file.path(basedir, "results")

segdir = file.path(progdir, "Segmentation")
source(file.path(segdir, "performance_functions.R"))

correct = "none"
# options = c("none", "N3", "N4", "N3_SS", "N4_SS",
#         "SyN", "SyN_sinc", "Rigid", "Affine", "Rigid_sinc", 
#         "Affine_sinc")
options = c("none", "N3_SS", "N4_SS", 
		"Rigid", "Rigid_sinc")

spec = matrix(c(
	'correct', 'c', 1, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

correct = opt$correct
print(opt)

# options = c("Rigid", "Rigid_sinc")

#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )

outfile = file.path(outdir, "111_Filenames.Rda")
load(file = outfile)


# iscen <- as.numeric(Sys.getenv("SGE_TASK_ID"))
# if (is.na(iscen)) iscen = 1


# scenarios = expand.grid(iimg = seq(nrow(fdf)), 
# 	correct = options, stringsAsFactors = FALSE )
# iimg = scenarios$iimg[iscen]
# correct = scenarios$correct[iscen]

iimg <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(iimg)) iimg = 71

runx = x = fdf[iimg,]


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

		# outfiles = nii.stub(basename(fdf$img))
		# outfiles = paste0(outfiles, "_predictors.Rda")
		# outfiles = file.path(fdf$outdir, outfiles)
		# file.exists(outfiles)
		# outfiles = nii.stub(basename(fdf$img))
		# outfiles = paste0(outfiles, "_models.Rda")
		# outfiles = file.path(fdf$outdir, outfiles)	
		# file.exists(outfiles)


	# fdf = fdf[1:10,]
	fpr.stop = .01
# run_model = function(x, fpr.stop = .1){
	fname= nii.stub(basename(x$img))
	fname = paste0(fname, "_predictors", adder, ".Rda")
	outfile = file.path(x$outdir, fname)
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
	mask.fname = switch(correct,
		"none" = x$mask,
		"N3" = x$mask,
		"N4" = x$mask,
		"N3_SS" = x$mask,
		"N4_SS" = x$mask,
		"SyN" = x$synssmask,
		"SyN_sinc" = x$sinc_synssmask,
		"Rigid" = x$rig_ssmask,
		"Affine" = x$aff_ssmask,
		"Rigid_sinc" = x$sinc_rig_ssmask,
		"Affine_sinc" = x$sinc_aff_ssmask
		)
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

	print(correct)

	xx = load(file=outfile)
	
	## Overlaid densities of ROIs
	# Rerun localization 
	df = img.pred$df
	# df$zscore_template = NULL
	df$mask = df$mask > 0
	stopifnot(all(df$Y %in% c(0, 1)))

	keep.ind = img.pred$keep.ind
	nim = img.pred$nim
	df = df[ keep.ind, ]
	miss.roi = img.pred$miss.roi
	rm(list="img.pred")

    keep.colnames = colnames(df)

	###########################################
	# Get cutoffs and make multiplier
	###########################################
	fname = file.path(outdir, 
        paste0("Aggregate_data_cutoffs", adder, ".Rda"))

    load(file = fname)
    
    df$gr_medztemp = TRUE
    keepnames = colnames(est.cutoffs)
    include = rep(TRUE, length=nrow(df))
    for (icut in keepnames){
        qcuts = est.cutoffs[, icut]
        colname = paste0(icut, ".cutoff")
        df[, colname] = df[, icut] >= qcuts[1] & 
            df[, icut] <= qcuts[2]
        include = include & df[, colname]
        print(icut)
    }



    sum(df$Y[!include])/ sum(df$Y)
    sum(include)/ length(include)
    df$include.all = include

    df$include = df$value >= 30 & df$value <= 100


    df$zval = df[, "zscore3.cutoff"] & df$include &
        df$pct_thresh.cutoff
    df$zval2 = df[, "zscore2.cutoff"] & df$zval
    df$zval_all = df[, "zscore_template.cutoff"] & 
        df$zval2

    ### Adding subject level above median
	med.ztemp = median(df$zscore_template)
	df$gr_medztemp = (df$zscore_template > med.ztemp)
	df$zval2_medztemp = df$zval2 &  df$gr_medztemp


    zval2 = nim
    zval2[keep.ind] = df$zval2
    zval2[is.na(zval2)] = 0
    zval2 = datatyper(zval2)
    zval2 = cal_img(zval2)
    outimg = file.path(x$outdir, 
    	paste0(nii.stub(fname, bn=TRUE), "_zval2"))
    writeNIfTI(zval2, filename= outimg)

    zval2 = nim
    zval2[keep.ind] = df$zval2_medztemp
    zval2[is.na(zval2)] = 0
    zval2 = datatyper(zval2)
    zval2 = cal_img(zval2)
    outimg = file.path(x$outdir, 
    	paste0(nii.stub(fname, bn=TRUE), "_zval2_medztemp"))
    writeNIfTI(zval2, filename= outimg)    

    # df$multiplier = df$zval2
    df$multiplier = df$zval2_medztemp
    df$candidate = df$multiplier | (df$Y == 1)
	#############################################
	# Case-Control Sampling
	#############################################
    seed = 20141022
    set.seed(seed)
    ich = which(df$Y == 1 & df$multiplier) 
    noich = which(df$Y != 1 & df$multiplier)

    mult.df = df[, !colnames(df) %in% keep.colnames]
    df = df[, colnames(df) %in% keep.colnames]
	# ich = which(df$Y == 1)
	# noich = which(df$Y != 1)

	size = 1e4
	prop = .25
	n.ich = ceiling(size*prop)
	n.noich = size - n.ich
	ich.ind = sample(ich, size=n.ich)
	noich.ind = sample(noich, size=n.noich)
	samp.ind = sort(c(ich.ind, noich.ind))

	# samp.ind = sample(nrow(df), size= 1e4)
	samps = seq(nrow(df)) %in% samp.ind
	train = df[samps,]

	fname= nii.stub(basename(x$img))
	fname = paste0(fname, "_predictors", adder, "_training.Rda")
	outfile = file.path(x$outdir, fname)

	mult.train = mult.df[samps,]

	save(train, mult.train, samps, keep.ind, file=outfile)

	test = df[!samps,]

	mult.test = mult.df[!samps,]

	pval = function(mod){
		cc =coef(summary(mod))[, 4]
		cc = cc[ !(names(cc) %in% "(Intercept)")]
	}

	runmod = function(formstr){
		form = as.formula(formstr)
		mod = glm(formula=form, data=train, family=binomial())
		return(mod)
	}
	formstr = "Y ~ . - mask - win_z"

	mod = runmod(formstr)
	takeout = colnames(train)
	takeout = takeout[ !(takeout %in% c("Y", "mask"))]
	take.mods = lapply( takeout , function(x) {
		runmod( formstr= paste0(formstr, " - ", x))
	})
	# pvals = sort(pval(mod.all))
	# top10 = names(pvals[1:10])

	test.pred = predict(mod, newdata=test, type="response")

	test.pred = test.pred * mult.test$multiplier
	smod = summary(mod)
	smod$deviance.resid = NULL
	smod$residuals = NULL
	mod = remove_lmparts(mod)



	######################################
	# Get all performance measures
	######################################
	pred <- prediction( test.pred, test$Y)
	perf <- performance(pred, "tpr", "fpr")
	pauc.cut = t(opt.cut(perf, pred))

	### Cutoff with maximal sensitivity constrained to fpr
    sens.cut = get_senscut(pred, fpr.stop=fpr.stop)

    #### Get cutoff with highest accuracy
    acc = get_acc(pred)
    print(acc)
    #### Get cutoff with highest pAUC
    pauc = get_pauc(pred, fpr.stop)
    print(pauc)
    #### Get cutoff with highest dice
	dice.coef = get_max_dice(pred)
	print(dice.coef)

	############################################
	# Prediction with taking out each individual model predictor
	############################################
	test.preds = sapply(take.mods, predict,
		newdata=test, type="response")
	iipred = 1

	test.preds = test.preds * mult.test$multiplier


	############################################
	# Creating Mtrices for output
	############################################
	lpred = length(take.mods)
	pauc.cuts = matrix(NA, nrow=lpred, ncol=ncol(pauc.cut))
	colnames(pauc.cuts) = colnames(pauc.cut)

	accs = matrix(NA, nrow=lpred, ncol=ncol(acc))
	colnames(accs) = colnames(acc)

	dice.coefs = accs
	colnames(dice.coefs) = colnames(dice.coef)

	sens.cuts = matrix(NA, nrow=lpred, ncol=ncol(sens.cut))
	colnames(sens.cuts) = colnames(sens.cut)

	paucs = rep(NA, length=lpred)

	pb = txtProgressBar(max=lpred, style=3)
	for (iipred in seq_along(take.mods)){
		tpred = test.preds[,iipred]

		ipred <- prediction( tpred, test$Y)
		# pred <- prediction( test.pred.all, test$Y)
		# pred <- prediction( test.pred.05, test$Y)
		iperf <- performance(ipred,"tpr","fpr")
		pauc.cuts[iipred,] = t(opt.cut(iperf, ipred))

		accs[iipred, ] = c(get_acc(ipred))

		paucs[iipred] = get_pauc(ipred, fpr.stop)

		dice.coefs[iipred,] = c(get_max_dice(ipred))

		sens.cuts[iipred,] = c(get_senscut(ipred, 
			fpr.stop=fpr.stop))

		setTxtProgressBar(pb, iipred)
	}
	close(pb)

	rownames(dice.coefs) = rownames(accs)= takeout
	rownames(sens.cuts) = rownames(dice.coefs)
	names(paucs) = takeout
	print(paucs)

	# Ys = matrix(test$Y, ncol=ncol(test.preds), 
	# 	nrow=nrow(test.preds))
	# preds <- prediction( test.preds, Ys)
	# perfs <- performance(preds,"tpr","fpr")
	# paucs.cut[ipred,] = t(opt.cut(perfs, preds))
	# paucs = performance(preds, "auc", fpr.stop= fpr.stop)


	# train.pred = predict(mod, newdata=train, type="response")


	# rownames(df) = NULL
	mods = list(mod=mod, 
		pauc = pauc, 
		train.ind = samps, 
		smod = smod)
	# return(mods)

	take.mods = lapply( take.mods , remove_lmparts)

	# mods = run_model(x)
	fname= nii.stub(basename(runx$img))
	fname = paste0(fname, "_models", adder, ".Rda")
	fname = file.path(x$outdir, fname)



    # gam.mod = gam(Y ~ 
    #     s(moment1) + 
    #     s(moment2) + 
    #     s(moment3) + 
    #     s(moment4) + 
    #     s(value) + 
    #     thresh +
    #     s(zscore1) + 
    #     s(zscore2) + 
    #     s(moment3) + 
    #     s(pct_thresh) + 
    #     pct_zero_neighbor + 
    #     any_zero_neighbor +
    #     s(dist_centroid) +
    #     s(smooth10) +
    #     s(smooth20)
    #     , data=train, family= binomial(), verbose=TRUE)

    # # gam.mod = remove_gamparts(gam.mod)

    # test.gam.pred = predict(gam.mod, test, type="response")
    # test.gam.pred = as.numeric(test.gam.pred)
    # cat("GAM Prediction \n")
    
    # gam.pred <- prediction( test.gam.pred, test$Y)
    # # pred <- prediction( test.pred.all, test$Y)
    # # pred <- prediction( test.pred.05, test$Y)
    # perf <- performance(gam.pred,"tpr","fpr")
    # gam.pauc.cut = t(opt.cut(perf, gam.pred))

    # gam.pauc = performance(gam.pred, "auc", fpr.stop= fpr.stop)
    # gam.pauc = gam.pauc@y.values[[1]] / fpr.stop
    # gam.pauc

    # gam.acc = performance(gam.pred, "acc")
    # ind = which.max(gam.acc@y.values[[1]])
    # gam.cutoff = gam.acc@x.values[[1]][[ind]]
    # gam.acc = gam.acc@y.values[[1]][ind]
    # if (is.infinite(gam.cutoff)){
    # 	gam.cutoff = 1
    # }
    # gam.acc = c(accuracy=gam.acc, cutoff= gam.cutoff)
    # gam.acc = t(gam.acc)
    # gam.acc


	save(mods, take.mods, 
		paucs, pauc, 
		sens.cut, sens.cuts,
		pauc.cut, pauc.cuts,
		dice.coef, dice.coefs,
		# gam.pauc, gam.acc, 
		# gam.pauc.cut, gam.mod,
		accs, acc, file = fname)


	# rm(list=c("perfs", "preds", "pred", 
	# 	"test.preds", "df", "perfs", "test", "tpred", 
	# 	"ipred", "iperf",
	# 	 "test.pred"))
	# for (i in 1:3) gc()

# }

# plot(perf)






