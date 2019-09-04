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
library(spm12r)
library(partykit)
library(magrittr)
set.seed(20150518)
homedir = "/Applications"
rootdir = file.path("/Volumes/DATA_LOCAL", 
    "Image_Processing")
if (Sys.info()[["user"]] %in% 
    "jmuschel") {
  homedir = "~"
  rootdir = file.path(
    "/legacy/dexter/disk2/smart", 
    "stroke_ct", "ident")
}
progdir = file.path(rootdir, 
    "programs")
basedir = file.path(rootdir, 
    "Registration")
tempdir = file.path(rootdir, 
    "Template")
atlasdir = file.path(tempdir, 
    "atlases")
outdir = file.path(basedir, 
    "results")

segdir = file.path(progdir, 
    "Reseg")
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
fdf$roi.fname = switch(correct,
    "none" = fdf$roi,
    "Rigid" = fdf$rig_ssroi
    )
fdf$img.fname = switch(correct,
        "none" = fdf$img,
        "Rigid" = fdf$rig_ssimg   
        )
fdf$res_rda_stub = file.path(fdf$outdir, 
    paste0("Reseg_Results_", 
        nii.stub(fdf$img.fname, bn = TRUE)
        )
    )

filename = file.path(outdir, 
    paste0("Reseg_Result_Formats", 
        adder, ".Rda"))
load(filename)

makedir = sapply( fdf$outdir, 
    function(x) {
    if (!file.exists(x)){
        dir.create(x, showWarnings =FALSE)
    }
})

##############################
# Keeping files where predictors exist
##############################
fdf$outfile = nii.stub(basename(fdf$img))
fdf$outfile = paste0("Reseg_", fdf$outfile, 
    "_predictors", adder, ".Rda")
fdf$outfile = file.path(fdf$outdir, 
    fdf$outfile)
fdf = fdf[file.exists(fdf$outfile), ]

###########################################
# Get cutoffs and make multiplier
###########################################
# mods = c("logistic", "lasso", 
#     "rf", 
#     "gam")
mods = "rf_reduced"
# , "gam", "cforest")
imod = mods[1]



iimg <- as.numeric(
    Sys.getenv("SGE_TASK_ID"))
if (is.na(iimg)) iimg = 21
rerun = TRUE

for (iimg in seq(nrow(fdf))) {
    x = fdf[iimg,]

    mask.fname = x$img.fname %>%
        nii.stub(bn = TRUE) %>%
        paste0(., 
            "_usemask.nii.gz") %>%
        file.path(x$preddir, .)

    usemask = readnii(mask.fname)
    img = readnii(x$img.fname)
    y = readnii(x$roi.fname)
    biny = y > 0.5

    mask = usemask > 0 | biny > 0
    vals = biny[mask]

    for (imod in mods){
        
        stub_fname = file.path(outdir, 
            paste0("Reseg_Aggregate_models", 
                adder))
        check_stub =  paste0(stub_fname, 
            "_", imod)
        check_fname = paste0(check_stub, 
            ".Rda")
        scheck_fname = paste0(check_stub, 
            "_smoothed.Rda")    

        res_rda = paste0(x$res_rda_stub, 
            "_", imod, ".Rda")

        if (
            (file.exists(check_fname) 
            & !file.exists(res_rda)) |
            rerun
            ) {

            load(check_fname)
            cutoff = modlist$mod.dice.coef[1,
                'cutoff']
            load(scheck_fname)
            scutoff = modlist$mod.dice.coef[1,
                'cutoff']            

            in_stub = x$img %>% 
                nii.stub(bn = TRUE) %>% 
                paste0("Reseg_", ., 
                    "_", imod) %>% 
                file.path(x$preddir, .)
            in_fname = paste0(in_stub, 
                    "_probability")
            pimg = readnii(in_fname)


            out_fname = paste0(in_stub, 
                    "_prediction")
            pred = datatyper(pimg > cutoff)
            writenii(pred, out_fname)

            out_fname = paste0(in_stub, 
                    "_prediction_cc")        
            cc = spm_bwlabel(pred, k = 100)
            writenii(cc, out_fname)


            sin_fname = paste0(in_fname, 
                "_smoothed")
            sm.pimg = readnii(sin_fname)

            out_fname = paste0(in_stub, 
                    "_smoothed_prediction")
            spred = datatyper(sm.pimg > 
                scutoff)
            writenii(spred, out_fname)

            out_fname = paste0(in_stub, 
                    "_smoothed_prediction_cc")  
            scc = spm_bwlabel(spred, k = 100)
            writenii(scc, out_fname)


            L = list(pred= pred,
                spred = spred, 
                cc = cc,
                scc = scc)
            results = lapply(L, 
                function(x){
                extrantsr:::sim(
                    vals,                    
                    x[mask])
            })
            save(results, file = res_rda)

        }
        print(imod)
    }


    print(iimg)
}