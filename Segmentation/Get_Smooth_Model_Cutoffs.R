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

    cut.filename = file.path(outdir, 
        paste0("Model_Cutoffs", adder, ".Rda"))
    x= load(cut.filename)


    fdf.run = fdf[run.ind, ]

    moddname = nii.stub(basename(fdf.run$img))
    moddname = file.path(fdf.run$outdir, 
        paste0(moddname, "_predictors", adder, ".Rda"))

    L = nrow(fdf.run)
    all.spreds = all.preds = all.df = NULL
    imod = 1

    for (imod in seq(L)){
        
        load(moddname[imod])

        df = img.pred$df
        keep.ind = img.pred$keep.ind
        nim = img.pred$nim

        id.outdir = fdf.run$outdir[imod]

        df$img = fdf.run$img[imod]

        #####################################
        # Load already-predicted images
        #####################################
        cn = names(all.cuts)
        icol = cn[1]
        preds = matrix(NA, nrow=nrow(df), 
            ncol=ncol(res)) 
        colnames(preds) = cn
        spreds = preds
        outstub = nii.stub(fdf.run$img[imod], bn=TRUE)
        for (icol in cn){
            outimg = file.path(id.outdir, 
                paste0(outstub, "_", icol, adder))
            img = readNIfTI(fname = outimg, reorient= FALSE)
            stopifnot(prod(dim(img)) == nrow(df))
            preds[, icol] = c(img)

            outimg = paste0(outimg, "_smoothed", adder)
            img = readNIfTI(fname = outimg, reorient= FALSE)
            stopifnot(prod(dim(img)) == nrow(df))
            spreds[, icol] = c(img)

            print(icol)
        }
        preds = preds[keep.ind,]
        spreds = spreds[keep.ind,]
        df = df[ keep.ind, ]

        ich = which(df$Y == 1)
        noich = which(df$Y != 1)
        df$mask = df$mask > 0

        size = 1e4
        prop = .25
        n.ich = ceiling(size*prop)
        n.noich = size - n.ich
        ich.ind = sample(ich, size=n.ich)
        noich.ind = sample(noich, size=n.noich)
        samp.ind = sort(c(ich.ind, noich.ind))

        df = df[samp.ind,]

        df$include = df$value >= 30 & df$value <= 100        
        preds = preds[samp.ind, ] * df$include 
        spreds = spreds[samp.ind, ] * df$include 

        all.df = rbind(all.df, df)
        all.spreds = rbind(all.spreds, spreds)
        all.preds = rbind(all.preds, preds)
        rm(list=c("img.pred", "df", "preds", "spreds"))
        print(imod)
    }



	#### only values 30 to 100 are likely to be ROI
	fpr.stop = 0.01

	opt.cut = function(perf, pred){
		cut.ind = mapply(FUN=function(x, y, p){
  			d = (x - 0)^2 + (y-1)^2
  			ind = which(d == min(d))
  			c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
  				cutoff = p[[ind]])
		}, perf@x.values, perf@y.values, pred@cutoffs)
	}


	r.sres = alply(all.spreds, 2, function(nd){
		pred <- prediction( nd, all.df$Y)				
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

	saccs = sapply(r.sres, `[[`, "acc")
	spaucs = sapply(r.sres, `[[`, "pauc")
	spauc.cuts = sapply(r.sres, `[[`, "pauc.cut")[3,]

	colnames(saccs) = names(spaucs) = names(spauc.cuts) = 
		colnames(all.spreds)

	all.scuts = saccs["cutoff", ]
	all.spauc.cuts = spauc.cuts

    names(all.spauc.cuts) = names(all.scuts) = cn

    all.scuts[is.infinite(all.scuts)] = 1
    all.spauc.cuts[is.infinite(all.spauc.cuts)] = 1

	mod.filename = file.path(outdir, 
		paste0("Smooth_Model_Cutoffs", adder, ".Rda"))

	save(all.scuts, all.spauc.cuts, file=mod.filename)

# }