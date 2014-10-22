###################################################################
## This code is for aggregate prediction of Image Segmentation of 
## CT
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

correct = "N3_SS"
options = c("none", "N3", "N4", "N3_SS", "N4_SS",
        "SyN", "SyN_sinc", "Rigid", "Affine")

short_predict = function(object, newdata, 
    lthresh=  .Machine$double.eps^0.5){
    tt <- terms(object)
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, 
        na.action = na.pass, 
        xlev = object$xlevels)
    if (is.null(cl <- attr(Terms, "dataClasses"))) 
            stop("no dataclasses")
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    # p <- object$rank
    beta <- object$coefficients
    beta = beta[ !is.na(beta) ]
    predictor = drop(X[, names(beta), drop=FALSE ] %*% beta)
   
    predictor <- family(object)$linkinv(predictor)
    predictor[ predictor < lthresh] = 0
    predictor
}

opt.cut = function(perf, pred){
    cut.ind = mapply(FUN=function(x, y, p){
        d = (x - 0)^2 + (y-1)^2
        ind = which(d == min(d))
        c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
            cutoff = p[[ind]])
    }, perf@x.values, perf@y.values, pred@cutoffs)
}



#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )

outfile = file.path(outdir, "111_Filenames.Rda")
load(file = outfile)

icorr <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(icorr)) icorr = 4
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
        "Affine" = "_Affine")

    filename = file.path(outdir, 
        paste0("Collapsed_Models", adder, ".Rda"))
    load(filename)

    makedir = sapply( fdf$outdir, function(x) {
    	if (!file.exists(x)){
    		dir.create(x, showWarnings =FALSE)
    	}
    })
    irow = 1
    x = fdf[irow,]

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



    all.df = NULL
    for (imod in seq(nrow(fdf.run))){
        load(moddname[imod])

        df = img.pred$df
        keep.ind = img.pred$keep.ind
        nim = img.pred$nim
        df = df[ keep.ind, ]

        df$img = fdf.run$img[imod]
        all.df = rbind(all.df, df)
        rm(list=c("img.pred", "df"))
        print(imod)
    }

    ich = which(all.df$Y == 1)
    noich = which(all.df$Y != 1)
    all.df$mask = all.df$mask > 0

    size = 1e5
    prop = .25
    n.ich = ceiling(size*prop)
    n.noich = size - n.ich
    ich.ind = sample(ich, size=n.ich)
    noich.ind = sample(noich, size=n.noich)
    samp.ind = sort(c(ich.ind, noich.ind))

    # samp.ind = sample(nrow(df), size= 1e4)
    samps = seq(nrow(all.df)) %in% samp.ind  

    # train = all.df[samps,]
    test = all.df[!samps,]
	test$include = test$value >= 30 & test$value <= 100

    fname = file.path(outdir, 
        paste0("Aggregate_data_cutoffs", adder, ".Rda"))

    load(file = fname)
    keepnames = colnames(est.cutoffs)
    include = rep(TRUE, length=ncol(df))
    for (icut in keepnames){
        qcuts = est.cutoffs[, icut]
        include = include & 
            (test[, icut] >= qcuts[1] & test[, icut] <= qcuts[2])
    }

    test$include.all = include

	preds = t(laply(all.mods, short_predict, newdata= test, 
                  .progress = "text"))

    rowMean = rowMeans(preds)
    rowMed = rowMedians(preds)
    rowMax = rowMaxs(preds)
    rowMin = rowMins(preds)
    rowProd = exp(rowSums(log(preds)))
    rowGeom = exp(rowMeans(log(preds)))

    preds = cbind(rowMean, rowMed, rowMax, rowMin,
        rowProd, rowGeom)
    rownames(preds) = NULL


	#### only values 0 to 100 are likely to be ROI
	# preds = preds * df$in0100
	preds = preds * test$include
	fpr.stop = 0.01

	opt.cut = function(perf, pred){
		cut.ind = mapply(FUN=function(x, y, p){
  			d = (x - 0)^2 + (y-1)^2
  			ind = which(d == min(d))
  			c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
  				cutoff = p[[ind]])
		}, perf@x.values, perf@y.values, pred@cutoffs)
	}


	r.res = alply(preds, 2, function(nd){
		pred <- prediction( nd, test$Y)				
		acc = performance(pred, "acc")
		ind = which.max(acc@y.values[[1]])
		cutoff = acc@x.values[[1]][ind]
		acc = acc@y.values[[1]][ind]
		acc = c(accuracy=acc, cutoff= cutoff)

		pauc = performance(pred, "auc", fpr.stop= fpr.stop)
		pauc = pauc@y.values[[1]] / fpr.stop
		
		perf <- performance(pred, "tpr", "fpr")
		pauc.cut = t(opt.cut(perf, pred))

		list(acc = acc, pauc = pauc, pauc.cut = pauc.cut)
	}, .progress = "text")

	accs = sapply(r.res, `[[`, "acc")
	paucs = sapply(r.res, `[[`, "pauc")
	pauc.cuts = sapply(r.res, `[[`, "pauc.cut")[3,]

	colnames(accs) = names(paucs) = names(pauc.cuts) = 
		colnames(preds)

	cuts = accs["cutoff", ]
	all.cuts = c(all.cutoffs, cuts, gam=gam.acc[, "cutoff"])
	all.pauc.cuts = c(all.pauc.cutoffs, pauc.cuts, 
        gam=gam.pauc.cut[, "cutoff"])
    names(all.pauc.cuts)[length(all.pauc.cuts)] = "gam"
    names(all.cuts)[length(all.cuts)] = "gam"

    all.cuts[is.infinite(all.cuts)] = 1
    all.pauc.cuts[is.infinite(all.pauc.cuts)] = 1

	mod.filename = file.path(outdir, 
		paste0("Model_Cutoffs", adder, ".Rda"))

	save(all.cuts, all.pauc.cuts, file=mod.filename)

# }