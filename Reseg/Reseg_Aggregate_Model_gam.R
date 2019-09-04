####################################
## This code is for aggregate models
## CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
####################################
#####################################
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
library(partykit)
library(parallel)
set.seed(20150518)
homedir = "/Applications"
rootdir = file.path("/Volumes/DATA_LOCAL", 
    "Image_Processing")
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = file.path("/dexter/disk2/smart", 
    "stroke_ct", "ident")
}
progdir = file.path(rootdir, "programs")
basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")
outdir = file.path(basedir, "results")

segdir = file.path(progdir, "Reseg")
source(file.path(segdir, 
    "Reseg_performance_functions.R"))

correct = "Rigid"

options = c("Rigid")
# options = c("none", "Rigid")


#### load voxel data
outfile = file.path(outdir, 
    "Reseg_111_Filenames.Rda")
load(file = outfile)

# for (correct in options){
    if ("all.df" %in% ls()){
        rm(list="all.df")
    }
    for (i in 1:3) gc()
    correct = match.arg(correct, options)
    adder = switch(correct, 
        "none"= "",
        "Rigid" = "_Rigid")


    filename = file.path(outdir, 
        paste0("Reseg_Result_Formats", 
            adder, ".Rda"))
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
    outfiles = paste0("Reseg_", outfiles, 
        "_predictors", adder, ".Rda")
    outfiles = file.path(fdf$outdir, 
        outfiles)
    fdf = fdf[file.exists(outfiles), ]

    
    fname = file.path(outdir, 
        paste0(
        "Reseg_Candidate_Aggregate_data", 
            adder, ".Rda"))
    load(fname)
    

    train$gr_medztemp = NULL
    train$mode = fdf$mode[match(train.img, 
        fdf$img)]
    test$mode = fdf$mode[match(test.img, 
        fdf$img)]

    runnames = colnames(train)
    nosmooth = c("any_zero_neighbor",
            "thresh", "pct_zero_neighbor")
    runnames = runnames[ !(runnames %in% 
        c("mask", "Y", "img", nosmooth))]


    runmod = function(formstr){
        form = as.formula(formstr)
        mod = glm(formula=form, 
            data=train, 
            family=binomial())
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

    test$multiplier = test.mult.df[, 
        "multiplier"]


    mod.time = system.time({
        mod = runmod(formstr)
    })
    mm = model.matrix(mod)
    mm = mm[ ,
        !colnames(mm) %in% c("(Intercept)")
        ]

    N = nrow(test)
    fpr.stop = 0.01

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
        control = list(trace = TRUE),
        method = "fREML")
    })
    # + mode
    gam.time


    ##############################
    # GAM PREDs
    ##############################
    test.gam.pred = rep(0, length=nrow(test))
    test.gam.pred[test$multiplier] = as.numeric(
        predict(gam.mod, 
            test[test$multiplier,], 
            type="response")
    )



    cat("GAM Prediction \n")
    
    test.gam.pred = as.numeric(test.gam.pred)
    test.gam.pred[ test.gam.pred > 1] = 1
    test.gam.pred = test.gam.pred * test$multiplier

    gam.modlist = mod_func(
        test.gam.pred, 
        test$Y, 
        fpr.stop)
    gam.modlist$mod = gam.mod
    gam.modlist$mod.time = gam.time
    gam.modlist$mod.perf = NULL

    modlist = gam.modlist

    stub_fname = file.path(outdir, 
        paste0("Reseg_Aggregate_models", 
            adder))
    fname = paste0(stub_fname, 
        "_gam.Rda")

    save(fdf.run, 
        modlist,  
        file = fname)

    

    # all.df$img = img

    # if (!is.null(cl)) stopCluster(cl)
