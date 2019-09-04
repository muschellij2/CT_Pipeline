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
library(MASS)
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
    

    train$gr_medztemp = NULL
    train$mode = fdf$mode[match(train.img, fdf$img)]
    test$mode = fdf$mode[match(test.img, fdf$img)]
    train$z = train$moment1 / train$moment2

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
    smod = summary(mod)
    mm = stepAIC(mod)

    mod = keep_mod(mod)
    test$z = test$moment1 / test$moment2
    test$multiplier = test.mult.df[, "multiplier"]

    test.pred = predict(mod, 
        test, 
        type = "response")
    # test.pred = test.pred * test$multiplier

    pred <- prediction( test.pred, test$Y)
    fpr.stop = .01
    # pred <- prediction( test.pred.all, test$Y)
    # pred <- prediction( test.pred.05, test$Y)
    perf <- performance(pred,"tpr","fpr")

    auc_pauc = function(pred, fpr.stop){
        auc = performance(pred, "auc")@y.values[[1]]
        pauc = performance(pred, "auc", 
            fpr.stop = fpr.stop)@y.values[[1]] / fpr.stop
        return(c(auc = auc, pauc = pauc))
    }    
    vals = auc_pauc(pred, fpr.stop = fpr.stop)
    print(vals)

    pred <- prediction( test.pred * test$multiplier, 
        test$Y)
    perf <- performance(pred,"tpr","fpr")   
    vals = auc_pauc(pred, fpr.stop = fpr.stop)
    print(vals)  
    # vals = auc_pauc(pred, fpr.stop = 0.05)

