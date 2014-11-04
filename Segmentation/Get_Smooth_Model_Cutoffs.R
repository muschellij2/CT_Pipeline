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

correct = "Rigid_sinc"
# options = c("none", "N3", "N4", "N3_SS", "N4_SS",
#         "SyN", "SyN_sinc", "Rigid", "Affine", "Rigid_sinc", 
#         "Affine_sinc")
options = c("none", "N3_SS", "N4_SS", 
      "Rigid", "Rigid_sinc")

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
if (is.na(icorr)) icorr = 5
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

    fname = file.path(outdir, 
        paste0("Aggregate_data", adder, ".Rda"))
    load(fname)

    all.df$mask = all.df$mask > 0
    all.df$multiplier = mult.df$multiplier

    all.ind = which(all.df$multiplier)

    all.preds = matrix(0, nrow=nrow(all.df), ncol= length(all.mods))
    for (imod in seq_along(all.mods)){
        x = short_predict(
            object = all.mods[[imod]],
            newdata=all.df[all.ind,])
        all.preds[all.ind, imod] = x 
        print(imod)
    }

    colnames(all.preds) = names(all.mods)
    # all.preds = t(laply(all.mods, short_predict, newdata= all.df, 
    #               .progress = "text"))
    add.preds = matrix(0, nrow=nrow(all.df), ncol= 7)
    colnames(add.preds) = c("rowMean", "rowMed", "rowMax", 
        "rowMin", "rowProd", "rowGeom", "gam")

    add.preds[ all.ind, "rowMean"] = 
        rowMeans(all.preds[  all.ind, ])
    add.preds[ all.ind, "rowMed"] = 
        rowMedians(all.preds[  all.ind, ])    
    add.preds[ all.ind, c("rowMin", "rowMax")] = 
        rowRanges(all.preds[  all.ind, ]) 

    nr = nrow(all.preds)
    add.preds[ all.ind, c("rowProd")] = 
        exp(rowSums(log(all.preds[all.ind, ])))
    add.preds[ all.ind, c("rowGeom")] = 
       add.preds[ all.ind, c("rowProd")] ^ (1/nr)


    gam.pred = predict(gam.mod, 
        all.df[all.ind,], 
        newdata.guaranteed = TRUE,
        type="response", block.size=5e5)
    cat("GAM Prediction \n")
    
    gam.pred = as.numeric(gam.pred)

    add.preds[all.ind, "gam"] = gam.pred

    all.preds = cbind(all.preds, add.preds)
    rm(list="add.preds")

    rownames(all.preds) = NULL

    all.preds = all.preds * all.df$multiplier

    # train = all.df[samps,]
    # test$include = test$value >= 30 & test$value <= 100

    L = nrow(fdf.run)
    all.spreds = array(NA, dim = dim(all.preds))
    imod = 1
    
    for (imod in seq(L)){
        
        imgname = fdf.run$img[imod]
        ind = which(img %in% imgname)
        keep.ind = l.keep.ind[[imod]]
        stopifnot(length(ind) == length(keep.ind))

        pred = all.preds[ind, ]

        load(moddname[imod])
        nim = img.pred$nim
        rm(list="img.pred")
        ipred = 1

        for (ipred in seq(ncol(all.spreds))) {
            orig.img = nim
            orig.img[keep.ind] = pred[, ipred]
            orig.img[is.na(orig.img)] = 0

            #### smooth the results
            sm.img  = mean_image(orig.img, nvoxels = 1)
            # sm.img = local_moment(img, nvoxels = 1, moment =1,
            #   center = FALSE)[[1]]
            sm.img = sm.img[keep.ind]
            sm.img[abs(sm.img) <  .Machine$double.eps^0.5 ] = 0
            all.spreds[ind, ipred] = sm.img
            # all.spreds[roi.not.in, ipred] = 0
            print(ipred)
        }
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

    cn = colnames(all.spreds)

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