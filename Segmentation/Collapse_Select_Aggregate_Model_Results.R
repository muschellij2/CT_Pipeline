###################################################################
## This code is for aggregate prediction of Image Segmentation of 
## CT, with all subsets
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
library(car)
library(getopt)
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

spec = matrix(c(
    'correct', 'c', 1, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

correct = opt$correct
print(opt)

options = c("none", "N3_SS", "N4_SS", 
      "Rigid", "Rigid_sinc")

if (is.null(correct)) correct = "Rigid"

for (correct in options){
    
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

   

    fname = file.path(outdir, 
            paste0("All_Subsets_Aggregate_Model", 
                adder, ".Rda")) 
    load(fname)

    all.scen = coef.mat != 0
    ################################
    # Making chunks for looping
    ################################
    chunksize = 300
    nchunks = ceiling(n.forms/chunksize)

    all.dice = all.acc = matrix(nrow=n.forms, ncol=2)
    all.pauc =  matrix(nrow=n.forms, ncol=1)
    all.senscut = matrix(nrow=n.forms, ncol=4)

    for (ichunk in seq(nchunks)){
        file = 
        file.path(outdir, 
            paste0("All_Subsets_Aggregate_Model_Results_Chunk_",
                ichunk, 
                adder, ".Rda"))
        xres = load(file)

        all.dice[ind, ]= dice
        all.acc[ind, ]= acc
        all.pauc[ind, ]= pauc
        all.senscut[ind, ]= senscut
        rm(list=xres)
        print(ichunk)
    }

    ichunk = 1
    file = 
    file.path(outdir, 
        paste0("All_Subsets_Aggregate_Model_Results_Chunk_",
            ichunk, 
            adder, ".Rda"))
    xres = load(file)

    colnames(all.dice) = colnames(dice)
    colnames(all.acc) = colnames(acc)
    colnames(all.senscut) = colnames(senscut)
    colnames(all.pauc) = "pauc"

    rm(list=xres)


    save(all.dice, all.acc, all.senscut, all.pauc,
        file = 
        file.path(outdir, 
            paste0("All_Subsets_Aggregate_Model_Results", 
                adder, ".Rda")))

    all.sensitivity = all.senscut
    all.dice = data.frame(all.dice)
    all.sensitivity = data.frame(all.sensitivity)
    all.acc = data.frame(all.acc)
    all.pauc = data.frame(all.pauc)

    measures = c("dice", "acc", "pauc", "sensitivity")
    all.meas = lapply(measures, 
        function(x){
        string = paste0("all.", x)
        obj = get(string)
        obj$value = obj[, x]
        obj = obj[, "value", drop=FALSE]
        ### -1 for intercept        
        obj$n_pred = rowSums(all.scen)-1
        obj$ind = seq(n.forms)  
        obj
    })
    names(all.meas) = measures

    bests = sapply(all.meas, function(x){
        which.max(x$value)
    })

    best.mods = all.scen[bests,]
    rownames(best.mods) = measures


    mapply(function(df, lab){
        plot(df$n_pred, df$value, ylab=lab, 
            xlab="Number of Predictors")
    }, all.meas, c("Dice", "Accuracy", "pAUC", "Sensitivity"))


    best.n.mods = lapply(all.meas, function(df){
        ddply(df, .(n_pred), function(x){
            ind = which.max(x$value)
            mod.ind = x$ind[ind]
            val = x$value[ind]
            c(model=mod.ind, value=val)
        })
    })


    sapply(best.n.mods, function(x){
        # x = x[x$n_pred > 0,]
        plot(value ~ n_pred, data=x)
    })


    var.imp = sapply(best.n.mods, function(x){
        colMeans(all.scen[x$model,])
    })

}