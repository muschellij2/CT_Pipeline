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

library(extrantsr)
library(randomForest)
library(methods)
library(glmnet)
set.seed(20150518)
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

correct = "Rigid"
# options = c("none", "N3", "N4", "N3_SS", "N4_SS",
#         "SyN", "SyN_sinc", "Rigid", "Affine", "Rigid_sinc", 
#         "Affine_sinc")
options = c("none", "N3_SS", "N4_SS", 
      "Rigid", "Rigid_sinc")

# options = c("none", "Rigid")


#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )


outfile = file.path(outdir, 
    "111_Filenames_with_volumes_stats.Rda")
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


    filename = file.path(outdir, 
        paste0("Result_Formats", adder, ".Rda"))
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

        paste0("Candidate_Aggregate_data", adder, ".Rda"))
    load(fname)
    
    # train$median = fdf$median[match(train.img, fdf$img)]
    # train$mean = fdf$mean[match(train.img, fdf$img)]
    train$gr_medztemp = NULL
    train$mode = fdf$mode[match(train.img, fdf$img)]
    test$mode = fdf$mode[match(test.img, fdf$img)]
    # train$rat = (train$value - train$mode)/train$mode * 1000

    runnames = colnames(train)
    nosmooth = c("any_zero_neighbor",
            "thresh", "pct_zero_neighbor")
    runnames = runnames[ !(runnames %in% 
        c("mask", "Y", "img", nosmooth))]



    runmod = function(formstr){
        form = as.formula(formstr)
        mod = glm(formula=form, data=train, family=binomial())
        return(mod)
    }

    sds = colSds(as.matrix(train))
    names(sds) = colnames(sds)
    novar = sds == 0
    novar = colnames(train[novar])
    novar = novar[!novar %in% c("mask")]
    novar = c("", novar)
    novar = paste0(novar, collapse=" - ")
    formstr = paste0("Y ~ . - mask", novar) 

    mod = runmod(formstr)
    mm = model.matrix(mod)
    mm = mm[ ,!colnames(mm) %in% c("(Intercept)")]

    takeout = colnames(train)
    takeout = takeout[ !(takeout %in% c("Y", "mask", "img"))]
    take.mods = llply( takeout , function(x) {
        runmod( formstr= paste0(formstr, " - ", x))
        }, .progress = "text")


    test$multiplier = test.mult.df[, "multiplier"]

    mod = remove_lmparts(mod)
    # gam.mod = remove_gamparts(gam.mod)
    take.mods = llply(take.mods, remove_lmparts, 
        .progress = "text")

    test.pred = short_predict(mod, test)
    test.pred = test.pred * test$multiplier

    pred <- prediction( test.pred, test$Y)
    fpr.stop = .01
    # pred <- prediction( test.pred.all, test$Y)
    # pred <- prediction( test.pred.05, test$Y)
    perf <- performance(pred,"tpr","fpr")
    pauc.cut = t(opt.cut(perf, pred))

    ###################################
    # Get measures and cutoffs
    ###################################
    N = nrow(test)
    sens.cut = get_senscut(pred, fpr.stop=fpr.stop) 
        # N = N,
        # predictions=test.pred, Y = test$Y)
        print(sens.cut)
    pauc = get_pauc(pred, fpr.stop=fpr.stop)
        print(pauc)
    dice.coef = get_max_dice(pred)
        print(dice.coef)
    acc = get_acc(pred)
        print(acc)

    ###################################
    # LASSO Model
    ###################################
    lasso.mod <- cv.glmnet(mm, train$Y, 
        type.measure = "class", 
        family = "binomial")

    test.mm = model.matrix(as.formula(formstr), test)
    test.mm = test.mm[, colnames(mm)]

    test.lasso.pred = predict(lasso.mod, 
        newx = test.mm,
        s = "lambda.1se", 
        type = "response")[, "1"]
    test.lasso.pred[!test$multiplier] = 0
    names(test.lasso.pred) = NULL
    test.lasso.pred[ test.lasso.pred > 1] = 1

    lasso.pred <- prediction( test.lasso.pred, test$Y)
    lasso.perf <- performance(lasso.pred,"tpr","fpr")
    lasso.pauc.cut = t(opt.cut(lasso.perf, lasso.pred))


    lasso.sens.cut = get_senscut(lasso.pred, 
        fpr.stop=fpr.stop)
        # N = N, predictions=test.lasso.pred, Y = test$Y)
        print(lasso.sens.cut)
    lasso.pauc = get_pauc(lasso.pred, fpr.stop=fpr.stop)
    print(lasso.pauc)
    lasso.dice.coef = get_max_dice(lasso.pred)
    print(lasso.dice.coef)
    lasso.acc = get_acc(lasso.pred)
    print(lasso.acc)    
    # lasso = glmnet(x = mm, y = train$Y, 
    #     family="binomial")


    ##############################
    # Random Forest Model
    ##############################
    rf.time <- system.time({
        rf.mod = randomForest(x = mm,
        y = factor(train$Y), 
        do.trace = TRUE)
        })

    test.rf.pred = rep(0, length=nrow(test))
    test.rf.pred[test$multiplier] = as.numeric(
        predict(rf.mod, test[test$multiplier,], 
            type="prob")[,"1"]
    )

    ##############################
    # Random Forest Predictions
    ##############################
    cat("# randomForest Prediction \n")
    
    test.rf.pred = as.numeric(test.rf.pred)
    test.rf.pred[ test.rf.pred > 1] = 1
    test.rf.pred = test.rf.pred * test$multiplier

    rf.pred <- prediction( test.rf.pred, test$Y)
    rf.perf <- performance(rf.pred,"tpr","fpr")
    rf.pauc.cut = t(opt.cut(rf.perf, rf.pred))


    rf.sens.cut = get_senscut(rf.pred, fpr.stop=fpr.stop)
        # N = N, predictions=test.rf.pred, Y = test$Y)
        print(rf.sens.cut)
    rf.pauc = get_pauc(rf.pred, fpr.stop=fpr.stop)
        print(rf.pauc)
    rf.dice.coef = get_max_dice(rf.pred)
        print(rf.dice.coef)
    rf.acc = get_acc(rf.pred)
        print(rf.acc)


    ### used gam - but bam supposedly faster
    gam.time = system.time({
        gam.mod = bam(Y ~ 
        s(moment1) + 
        s(moment2) + 

        s(skew) + 
        s(kurtosis) + 
        s(value) + 
        thresh +
        s(zscore1) + 
        s(zscore2) + 

        s(zscore3) + 
        s(pct_thresh) + 
        pct_zero_neighbor + 
        any_zero_neighbor +
        s(dist_centroid) +
        s(smooth10) +

        s(smooth20) 
        , data=train, family= binomial(), 
        method = "fREML")
    })
    # + mode
    gam.time

    ##############################
    # GAM PREDs
    ##############################
    test.gam.pred = rep(0, length=nrow(test))
    test.gam.pred[test$multiplier] = as.numeric(
        predict(gam.mod, test[test$multiplier,], type="response")
    )

    cat("GAM Prediction \n")
    
    test.gam.pred = as.numeric(test.gam.pred)
    test.gam.pred[ test.gam.pred > 1] = 1
    test.gam.pred = test.gam.pred * test$multiplier

    gam.pred <- prediction( test.gam.pred, test$Y)

    perf <- performance(gam.pred,"tpr","fpr")
    gam.pauc.cut = t(opt.cut(perf, gam.pred))


    gam.sens.cut = get_senscut(gam.pred, fpr.stop=fpr.stop)
        # N = N, predictions=test.gam.pred, Y = test$Y)
        print(gam.sens.cut)
    gam.pauc = get_pauc(gam.pred, fpr.stop=fpr.stop)
        print(gam.pauc)
    gam.dice.coef = get_max_dice(gam.pred)
        print(gam.dice.coef)
    gam.acc = get_acc(gam.pred)
        print(gam.acc)


    iipred = 1
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
        imod = take.mods[[iipred]]
        tpred = short_predict(imod, newdata=test)
        tpred = tpred * test$multiplier


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
            # N = length(tpred),
            # predictions = tpred,
            # Y = test$Y))

        setTxtProgressBar(pb, iipred)   
    }
    close(pb)

    rownames(dice.coefs) = rownames(accs)= takeout
    rownames(sens.cuts) = rownames(dice.coefs)
    names(paucs) = takeout
    print(paucs)

    # all.df$img = img

    # rownames(df) = NULL
    mods = list(mod=mod, 

        pauc = pauc)
        # return(mods)
    # }


    fname = file.path(outdir, 
        paste0("Aggregate_models", adder, ".Rda"))

    save(mods, take.mods, paucs, pauc, fdf.run, 

        pauc.cut, pauc.cuts,
        sens.cut, sens.cuts,
        dice.coef, dice.coefs,
        accs, acc, 
        gam.sens.cut, 
        gam.time,
        gam.dice.coef,
        gam.mod, gam.acc, gam.pauc, 
        gam.pauc.cut,
        rf.time,
        rf.sens.cut, 
        rf.dice.coef,
        rf.mod, rf.acc, rf.pauc, 
        rf.pauc.cut,

        lasso.sens.cut, 
        lasso.dice.coef,
        lasso.mod, lasso.acc, lasso.pauc, 
        lasso.pauc.cut,        
        file = fname)
# }
