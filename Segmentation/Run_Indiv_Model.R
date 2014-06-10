#####################################################################
## This code is for Image Segmentation of CT
## The code is R but calls out FSL
##
## Author: John Muschelli
## Last updated: May 20, 2014
#####################################################################
#####################################################################
rm(list=ls())
library(plyr)
library(cttools)
library(devtools)
library(ROCR)
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

rerun = FALSE
outdir = file.path(basedir, "results")


#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )

fnames = gsub("^bws", "", fnames)
fnames = paste0(fnames, ".gz")
ids = gsub("(\\d\\d\\d-(\\d|)\\d\\d\\d)_.*", "\\1", fnames)
fdf = data.frame(id = ids, stringsAsFactors= FALSE)
fdf$iddir = file.path(basedir, fdf$id)
fdf$outdir = file.path(fdf$iddir, "Predictors")
makedir = sapply( fdf$outdir, dir.create, showWarnings =FALSE)
fdf$roi = file.path(rootdir, "ROI_data", fdf$id, fnames)
fdf$img = file.path(fdf$iddir, gsub("ROI\\.nii", ".nii", fnames))
fdf$mask = file.path(fdf$iddir, 
	"Skull_Stripped", 
	gsub("ROI\\.nii", "_SS_Mask_0.01.nii", fnames))
iimg <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(iimg)) iimg = 2
x = fdf[iimg,]

	# outfiles = nii.stub(basename(fdf$img))
	# outfiles = paste0(outfiles, "_predictors.Rda")
	# outfiles = file.path(fdf$outdir, outfiles)
	# file.exists(outfiles)
	# outfiles = nii.stub(basename(fdf$img))
	# outfiles = paste0(outfiles, "_models.Rda")
	# outfiles = file.path(fdf$outdir, outfiles)	
	# file.exists(outfiles)

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



# fdf = fdf[1:10,]
fpr.stop = .1
# run_model = function(x, fpr.stop = .1){

	fname= nii.stub(basename(x$img))
	fname = paste0(fname, "_predictors.Rda")
	outfile = file.path(x$outdir, fname)
	if (!file.exists(outfile) | rerun){
	system.time({
		img.pred = make_predictors(
		img=x$img, 
		mask = x$mask, 
		roi = x$roi,
		nvoxels = 1, 
		moments = 1:4, lthresh = 40, uthresh = 80,
		save_imgs = TRUE, 
		outdir = x$outdir,
		overwrite = FALSE, 
		verbose= TRUE)
	})
		save(img.pred, file=outfile, compress = TRUE)
	} else {
		load(file=outfile)
	}

	df = img.pred$df
	keep.ind = img.pred$keep.ind
	nim = img.pred$nim
	df = df[ keep.ind, ]
	miss.roi = img.pred$miss.roi


	# df$Y = roi[keep.ind]

	# img = readNIfTI(x$img, reorient= FALSE)
	# # mask = readNIfTI(x$mask, reorient= FALSE)
	# # erode the mask
	# mask = fslerode(file=x$mask, kopts = "-kernel box 1x1x1", 
	#                   reorient=FALSE, retimg = TRUE)
	# roi = readNIfTI(x$roi, reorient = FALSE)
	# mask = mask > 0
	# mask[ roi == 1] = 1

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
	formstr = "Y ~ ."

	mod = runmod(formstr)
	takeout = colnames(train)
	takeout = takeout[ !(takeout %in% "Y")]
	take.mods = lapply( takeout , function(x) {
		runmod( formstr= paste0(formstr, " - ", x))
	})
	# pvals = sort(pval(mod.all))
	# top10 = names(pvals[1:10])

	# mod.10 = glm(Y ~ ., data=train[, c("Y", top10)], 
	# 	family=binomial())

	# pvals = sort(pval(mod.all))
	# l.05 = names(which(pvals < .05))

	# mod.05 = glm(Y ~ ., data=train[, c("Y", l.05)], 
	# 	family=binomial())	



	test.pred = predict(mod, newdata=test, type="response")

	# train.pred = predict(mod, newdata=train, type="response")

	pred <- prediction( test.pred, test$Y)
	fdr.stop = .1
	# pred <- prediction( test.pred.all, test$Y)
	# pred <- prediction( test.pred.05, test$Y)
	perf <- performance(pred,"tpr","fpr")
	xind = perf@x.values[[1]] <= fdr.stop
	perf@x.values[[1]] = perf@x.values[[1]][xind]
	perf@y.values[[1]] = perf@y.values[[1]][xind]

	# auc = performance(pred, "auc")@y.values[[1]]
	# plot(perf)
	pauc = performance(pred, "auc", fpr.stop= fpr.stop)
	pauc = pauc@y.values[[1]] / fpr.stop
	pauc



	test.preds = sapply(take.mods, predict,
		newdata=test, type="response")

	Ys = matrix(test$Y, ncol=ncol(test.preds), nrow=nrow(test.preds))

	# train.pred = predict(mod, newdata=train, type="response")

	preds <- prediction( test.preds, Ys)
	# pred <- prediction( test.pred.all, test$Y)
	# pred <- prediction( test.pred.05, test$Y)
	# perf <- performance(pred,"tpr","fpr")
	# xind = perf@x.values[[1]] <= fdr.stop
	# perf@x.values[[1]] = perf@x.values[[1]][xind]
	# perf@y.values[[1]] = perf@y.values[[1]][xind]

	# auc = performance(pred, "auc")@y.values[[1]]
	# plot(perf)
	paucs = performance(preds, "auc", fpr.stop= fpr.stop)
	paucs = unlist(paucs@y.values) / fpr.stop
	paucs


	smod = summary(mod)
	smod$deviance.resid = NULL
	mod = remove_lmparts(mod)

	# rownames(df) = NULL
	mods = list(mod=mod, 
		pauc = pauc, 
		train.ind = samps, 
		smod = smod)
	# return(mods)
# }

# mods = run_model(x)
fname= nii.stub(basename(x$img))
fname = paste0(fname, "_models.Rda")
fname = file.path(x$outdir, fname)

save(mods, take.mods, paucs, pauc, file = fname)


# plot(perf)







# # function(x){
# irow = 5
# x = fdf[irow,]
# new.img = readNIfTI(x$img, reorient= FALSE)
# new.mask = readNIfTI(x$mask, reorient= FALSE)
# # erode the mask
# new.mask = fslerode(file=new.mask, kopts = "-kernel box 1x1x1", 
#                   reorient=FALSE, retimg = TRUE)

# new.roi = readNIfTI(x$roi, reorient = FALSE)

# ind = which(new.roi > 0 & new.mask < 1, arr.ind = TRUE)

# new.img.pred = make_predictors(img=new.img, mask = new.mask, nvoxels = 1, moments = 1:4,
#                             lthresh = 40, uthresh = 80)
# new.df = new.img.pred
# new.df = data.frame(new.df)
# new.df$Y = c(new.roi)

# new.pred = predict(mod, newdata=new.df, type="response")
# pimg = new.roi
# pimg@.Data = array(new.pred, dim=dim(pimg))
# pimg = cal_img(pimg)
# pimg[new.mask < 1 ] = 0

# run.pred = new.pred[new.mask > 0 | new.df$Y > 0]
# Y = new.df$Y[new.mask > 0 | new.df$Y > 0]
# # centroid = which(mask > 0, arr.ind=TRUE)
# # dcent = t(t(centroid) - colMeans(centroid))
# # dcent = sqrt(rowSums(dcent^2))

# pred <- prediction( run.pred, Y)
# perf <- performance(pred,"tpr","fpr")


# auc = performance(pred, "auc")@y.values[[1]]
# # plot(perf)
# fpr.stop = .1
# pauc = performance(pred, "auc", fpr.stop= fpr.stop)
# pauc = pauc@y.values[[1]] / fpr.stop
# pauc



