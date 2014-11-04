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
# 		"Rigid", "Rigid_sinc")
options = c("Rigid", "Rigid_sinc")

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
if (is.na(iimg)) iimg = 1

runx = x = fdf[iimg,]

remove_lmparts = function(mod){
	keep = c("coefficients", "xlevels", 
		"contrasts", "family", "terms")
	mn = names(mod)
	mn = mn[ !( mn %in% keep)]
	for (iname in mn){
		mod[[iname]] = NULL
	}
	mod
}

sub_samp = function(x, pct = .1, maxcut = 5e3){
	rr = range(x)
	breaks = seq(rr[1], rr[2], length.out = 10)
	cuts = cut(x, breaks= breaks, include.lowest=TRUE)
	levs = levels(cuts)
	l = lapply(levs, function(ilev){
		which(cuts == ilev)
	})
	n.inlev = sapply(l, length)
	n.inlev = ceiling(n.inlev*pct)
	n.inlev = pmin(n.inlev, maxcut)
	ind = sort(unlist(mapply(function(ind, n){
		sample(ind, size = n)
	}, l, n.inlev)))
	return(ind)
}


remove_gamparts = function(mod){
    keep = c("assign", "cmX", "coefficients", "contrasts", 
        "family", 
        "formula", "model", "na.action", "nsdf", "pred.formula", 
        "pterms", "smooth",  "Xcentre", "xlevels", "terms")
    # "Vp",
    mn = names(mod)
    mn = mn[ !( mn %in% keep)]
    for (iname in mn){
        mod[[iname]] = NULL
    }
    mod$model = mod$model[1,]
    mod
}



for (correct in options){
    
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
	fpr.stop = 	fdr.stop = .01
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
	df$mask = df$mask > 0
	stopifnot(all(df$Y %in% c(0, 1)))

	keep.ind = img.pred$keep.ind
	nim = img.pred$nim
	df = df[ keep.ind, ]
	miss.roi = img.pred$miss.roi
	rm(list="img.pred")

	#############################################
	# Case-Control Sampling
	#############################################
	ich = which(df$Y == 1)
	noich = which(df$Y != 1)

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
	test = df[!samps,]

	pval = function(mod){
		cc =coef(summary(mod))[, 4]
		cc = cc[ !(names(cc) %in% "(Intercept)")]
	}

	runmod = function(formstr){
		form = as.formula(formstr)
		mod = glm(formula=form, data=train, family=binomial())
		return(mod)
	}
	formstr = "Y ~ . - mask"

	mod = runmod(formstr)
	takeout = colnames(train)
	takeout = takeout[ !(takeout %in% c("Y", "mask"))]
	take.mods = lapply( takeout , function(x) {
		runmod( formstr= paste0(formstr, " - ", x))
	})
	# pvals = sort(pval(mod.all))
	# top10 = names(pvals[1:10])

	test.pred = predict(mod, newdata=test, type="response")

	smod = summary(mod)
	smod$deviance.resid = NULL
	smod$residuals = NULL
	mod = remove_lmparts(mod)	

	# train.pred = predict(mod, newdata=train, type="response")
	opt.cut = function(perf, pred){
		cut.ind = mapply(FUN=function(x, y, p){
  			d = (x - 0)^2 + (y-1)^2
  			ind = which(d == min(d))
  			c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
  				cutoff = p[[ind]])
		}, perf@x.values, perf@y.values, pred@cutoffs)
	}


	pred <- prediction( test.pred, test$Y)
	# pred <- prediction( test.pred.all, test$Y)
	# pred <- prediction( test.pred.05, test$Y)
	perf <- performance(pred, "tpr", "fpr")
	pauc.cut = t(opt.cut(perf, pred))
	xind = perf@x.values[[1]] <= fdr.stop
	perf@x.values[[1]] = perf@x.values[[1]][xind]
	perf@y.values[[1]] = perf@y.values[[1]][xind]

	# auc = performance(pred, "auc")@y.values[[1]]
	# plot(perf)
	pauc = performance(pred, "auc", fpr.stop= fpr.stop)
	pauc = pauc@y.values[[1]] / fpr.stop
	pauc

	acc = performance(pred, "acc")
	ind = which.max(acc@y.values[[1]])
	cutoff = acc@x.values[[1]][[ind]]
	acc = acc@y.values[[1]][ind]
	acc = t(as.matrix(c(accuracy=acc, cutoff= cutoff)))
	acc


	test.preds = sapply(take.mods, predict,
		newdata=test, type="response")
	iipred = 1
	lpred = length(take.mods)
	paucs.cut = matrix(NA, nrow=lpred, ncol=3)
	colnames(paucs.cut) = colnames(pauc.cut)

	accs = matrix(NA, nrow=lpred, ncol=2)
	colnames(accs) = colnames(acc)

	paucs = rep(NA, length=lpred)
	pb = txtProgressBar(max=lpred, style=3)
	for (iipred in seq_along(take.mods)){
		tpred = test.preds[,iipred]

		ipred <- prediction( tpred, test$Y)
		# pred <- prediction( test.pred.all, test$Y)
		# pred <- prediction( test.pred.05, test$Y)
		iperf <- performance(ipred,"tpr","fpr")
		paucs.cut[iipred,] = t(opt.cut(iperf, ipred))

		# auc = performance(pred, "auc")@y.values[[1]]
		# plot(perf)
		run.pauc = performance(ipred, "auc", 
			fpr.stop= fpr.stop)
		paucs[iipred] = unlist(run.pauc@y.values) / fpr.stop

		run.acc = performance(ipred, "acc")
		ind = sapply(run.acc@y.values, which.max)
		cutoffs = mapply(function(xval, ind){
			xval[[ind]]
		}, run.acc@x.values, ind)
		run.accs = mapply(function(xval, ind){
			xval[[ind]]
		}, run.acc@y.values, ind)
		accs[iipred, ] = c(accuracy=run.accs, cutoff= cutoffs)

		setTxtProgressBar(pb, iipred)	
	}
	close(pb)	

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
		pauc.cut, paucs.cut,
		# gam.pauc, gam.acc, 
		# gam.pauc.cut, gam.mod,
		accs, acc, file = fname)


	# rm(list=c("perfs", "preds", "pred", 
	# 	"test.preds", "df", "perfs", "test", "tpred", 
	# 	"ipred", "iperf",
	# 	 "test.pred"))
	# for (i in 1:3) gc()

}

# plot(perf)






