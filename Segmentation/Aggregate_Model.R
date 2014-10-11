#####################################################################
## This code is for prediciton of Image Segmentation of CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
#####################################################################
#####################################################################
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

correct = "SyN"
options = c("none", "N3", "N4", "N3_SS", "N4_SS",
        "SyN", "SyN_sinc", "Rigid", "Affine")
for (correct in options){
    rm(list="all.df")
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

    #### load voxel data
    outfile = file.path(outdir, "Voxel_Info.Rda")
    load(file=outfile )

    outfile = file.path(outdir, "111_Filenames.Rda")
    load(file = outfile)


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
    lmod = 10
    fdf.run = fdf[seq(lmod), ]

    moddname = nii.stub(basename(fdf.run$img))
    moddname = file.path(fdf.run$outdir, 
        paste0(moddname, "_predictors", adder, ".Rda"))

    all.df = NULL
    for (imod in seq(lmod)){
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

    size = 1e5
    prop = .25
    n.ich = ceiling(size*prop)
    n.noich = size - n.ich
    ich.ind = sample(ich, size=n.ich)
    noich.ind = sample(noich, size=n.noich)
    samp.ind = sort(c(ich.ind, noich.ind))

    # samp.ind = sample(nrow(df), size= 1e4)
    samps = seq(nrow(all.df)) %in% samp.ind
    
    img = all.df$img
    all.df$img = NULL    

    train = all.df[samps,]
    test = all.df[!samps,]


    runmod = function(formstr){
        form = as.formula(formstr)
        mod = glm(formula=form, data=train, family=binomial())
        return(mod)
    }
    formstr = "Y ~ . - mask"

    mod = runmod(formstr)
    takeout = colnames(train)
    takeout = takeout[ !(takeout %in% c("Y", "mask", "img"))]
    take.mods = llply( takeout , function(x) {
        runmod( formstr= paste0(formstr, " - ", x))
        }, .progress = "text")

    # gam.mod = gam(Y ~ 
    #     s(moment1) + 
    #     s(moment2) + 
    #     s(moment3) + 
    #     s(moment4) + 
    #     s(value) + 
    #     thresh
    #     s(zscore1) + 
    #     s(zscore2) + 
    #     s(moment3) + 
    #     s(pct_thresh) + 
    #     pct_zero_neighbor + 
    #     any_zero_neighbor +
    #     s(dist_centroid) +
    #     s(smooth10) +
    #     s(smooth20)
    #     , data=all.df)

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

    mod = remove_lmparts(mod)

    take.mods = llply(take.mods, remove_lmparts, .progress = "text")

    test.pred = short_predict(mod, test)

    pred <- prediction( test.pred, test$Y)
    fpr.stop = fdr.stop = .1
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



    # test.preds = laply(take.mods, short_predict,
    #     newdata=test, .progress = "text")
    i = 1
    paucs = rep(NA, length = length(take.mods))
    for (i in seq_along(take.mods)){
        imod = take.mods[[i]]
        test.preds = short_predict(imod, newdata=test)
        preds <- prediction( test.preds, test$Y)
        myperf = performance(preds, "auc", fpr.stop = fpr.stop)
        myperf = unlist(myperf@y.values) / fpr.stop
        paucs[i] = myperf
        print(i)
    }
    # Ys = matrix(test$Y, ncol=ncol(test.preds), nrow=nrow(test.preds))

    # train.pred = predict(mod, newdata=train, type="response")

    # preds <- prediction( test.preds, Ys)


        # auc = performance(pred, "auc")@y.values[[1]]
        # plot(perf)
    # paucs = performance(preds, "auc", fpr.stop= fpr.stop)
    # paucs = unlist(paucs@y.values) / fpr.stop
    # paucs

    # all.df$img = img

    # rownames(df) = NULL
    mods = list(mod=mod, 
        pauc = pauc, 
        train.ind = samps)
        # return(mods)
    # }


    fname = file.path(outdir, 
        paste0("Aggregate_models", adder, ".Rda"))

    save(mods, take.mods, paucs, pauc, fdf.run, 
        file = fname)
}
