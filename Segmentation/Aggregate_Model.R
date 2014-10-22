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

correct = "SyN"
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
    # fdf.run = fdf[run.ind, ]

    # moddname = nii.stub(basename(fdf.run$img))
    # moddname = file.path(fdf.run$outdir, 
    #     paste0(moddname, "_predictors", adder, ".Rda"))

    # all.df = NULL
    # for (imod in seq(nrow(fdf.run))){
    #     load(moddname[imod])

    #     df = img.pred$df
    #     keep.ind = img.pred$keep.ind
    #     nim = img.pred$nim
    #     df = df[ keep.ind, ]

    #     df$img = fdf.run$img[imod]
    #     all.df = rbind(all.df, df)
    #     rm(list=c("img.pred", "df"))
    #     print(imod)
    # }

    fname = file.path(outdir, 
        paste0("Aggregate_data", adder, ".Rda"))
    load(fname)
    
    runnames = colnames(all.df)
    nosmooth = c("any_zero_neighbor",
            "thresh", "pct_zero_neighbor")
    runnames = runnames[ !(runnames %in% 
        c("mask", "Y", "img", nosmooth))]

    fname = file.path(outdir, 
        paste0("Aggregate_data_cutoffs", adder, ".Rda"))

    load(file = fname)
    
    train = all.df[samps,]
    test = all.df[!samps,]



    keepnames = colnames(est.cutoffs)
    include = rep(TRUE, length=nrow(all.df))
    for (icut in keepnames){
        qcuts = est.cutoffs[, icut]
        include = include & 
            (all.df[, icut] >= qcuts[1] & 
                all.df[, icut] <= qcuts[2])
        print(icut)
    }

    sum(all.df$Y[!include])/ sum(all.df$Y)
    sum(include)/ length(include)
    all.df$include.all = include

    test$include.all = include[!samps]

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


    ### used gam - but bam supposedly faster
    gam.time = system.time({
        gam.mod = bam(Y ~ 
        s(moment1) + 
        s(moment2) + 
        s(moment3) + 
        s(moment4) + 
        s(value) + 
        thresh +
        s(zscore1) + 
        s(zscore2) + 
        s(moment3) + 
        s(pct_thresh) + 
        pct_zero_neighbor + 
        any_zero_neighbor +
        s(dist_centroid) +
        s(smooth10) +
        s(smooth20)
        , data=train, family= binomial(), 
        method = "fREML")
    })
    gam.time

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

    remove_gamparts = function(mod){
        keep = c("assign", "cmX", "coefficients", "contrasts", 
            "family", 
            "formula", "model", "na.action", "nsdf", 
            "pred.formula", 
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

    mod = remove_lmparts(mod)
    # gam.mod = remove_gamparts(gam.mod)
    take.mods = llply(take.mods, remove_lmparts, 
        .progress = "text")

    test.pred = short_predict(mod, test)
    test.pred = test.pred * test$include.all

    pred <- prediction( test.pred, test$Y)
    fpr.stop = fdr.stop = .01
    # pred <- prediction( test.pred.all, test$Y)
    # pred <- prediction( test.pred.05, test$Y)
    perf <- performance(pred,"tpr","fpr")
    pauc.cut = t(opt.cut(perf, pred))
    xind = perf@x.values[[1]] <= fdr.stop
    perf@x.values[[1]] = perf@x.values[[1]][xind]
    perf@y.values[[1]] = perf@y.values[[1]][xind]

    # auc = performance(pred, "auc")@y.values[[1]]
    # plot(perf)
    pauc = performance(pred, "auc", fpr.stop= fpr.stop)
    pauc = pauc@y.values[[1]] / fpr.stop
    pauc


    ##############################
    # GAM PREDs
    ##############################
    test.gam.pred = predict(gam.mod, test, type="response")
    cat("GAM Prediciton \n")
    
    test.gam.pred = as.numeric(test.gam.pred)
    test.gam.pred = test.gam.pred * test$include.all

    gam.pred <- prediction( test.gam.pred, test$Y)
    # pred <- prediction( test.pred.all, test$Y)
    # pred <- prediction( test.pred.05, test$Y)
    perf <- performance(gam.pred,"tpr","fpr")
    gam.pauc.cut = t(opt.cut(perf, gam.pred))

    gam.pauc = performance(gam.pred, "auc", fpr.stop= fpr.stop)
    gam.pauc = gam.pauc@y.values[[1]] / fpr.stop
    gam.pauc

    gam.acc = performance(gam.pred, "acc")
    ind = which.max(gam.acc@y.values[[1]])
    gam.cutoff = gam.acc@x.values[[1]][[ind]]
    gam.acc = gam.acc@y.values[[1]][ind]
    gam.acc = c(accuracy=gam.acc, cutoff= gam.cutoff)
    gam.acc = t(gam.acc)
    gam.acc



    # test.preds = laply(take.mods, short_predict,
    #     newdata=test, .progress = "text")
    i = 1
    lpred = length(take.mods)
    paucs = rep(NA, length = lpred)
    accs = matrix(NA, nrow=lpred, ncol=2)
    colnames(accs) = c("accuracy", "cutoff")

    paucs.cut = matrix(NA, nrow=lpred, ncol=3)
    colnames(paucs.cut) = colnames(pauc.cut)    
    for (i in seq_along(take.mods)){
        imod = take.mods[[i]]
        test.preds = short_predict(imod, newdata=test)
        test.preds = test.preds * test$include.all

        preds <- prediction( test.preds, test$Y)
        myperf = performance(preds, "auc", fpr.stop = fpr.stop)
        myperf = unlist(myperf@y.values) / fpr.stop
        paucs[i] = myperf

        #### get cutoff with highest accuracy
        acc = performance(preds, "acc")
        ind = which.max(acc@y.values[[1]])
        cutoff = acc@x.values[[1]][[ind]]
        acc = acc@y.values[[1]][ind]
        acc = c(accuracy=acc, cutoff= cutoff)
        accs[i, ] = acc

        iperf <- performance(preds,"tpr","fpr")
        paucs.cut[i,] = t(opt.cut(iperf, preds))

        print(i)
    }

    acc = performance(pred, "acc")
    ind = which.max(acc@y.values[[1]])
    cutoff = acc@x.values[[1]][[ind]]
    acc = acc@y.values[[1]][ind]
    acc = c(accuracy=acc, cutoff= cutoff)
    acc = t(acc)
    acc 

    

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
        pauc.cut, paucs.cut,
        gam.time,
        accs, acc, gam.mod, gam.acc,  gam.pauc, 
        gam.pauc.cut,
        file = fname)
# }
