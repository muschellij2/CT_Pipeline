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

mods = c("logistic", "lasso", "rf")
# , "gam", "cforest")
imod = mods[3]

fname = file.path(outdir, 
        paste0("Reseg_Aggregate_data", 
            adder, ".Rda"))
load(fname)

all.df$mode = fdf$mode[match(img, 
    fdf$img)]
all.df$mask = all.df$mask > 0
all.df$multiplier = mult.df$multiplier


fdf.run = fdf.run[ 
    fdf.run$group %in% "Train",]
L = nrow(fdf.run)
coll.mods = seq(L)

fdf.run$modname = nii.stub(fdf.run$img, 
    bn=TRUE)
fdf.run$modname = file.path(fdf.run$outdir, 
        paste0(fdf.run$modname, 
            "_predictors", 
            adder, ".Rda"))

# mods = c("logistic", "lasso", 
#     "rf", 
#     "gam",
#     "rf_reduced")
mods = "rf_reduced"
# , "gam", "cforest")
imod = mods[1]


stub_fname = file.path(outdir, 
    paste0("Reseg_Aggregate_models", 
        adder))

for (imod in mods){
    check_stub =  paste0(stub_fname, 
        "_", imod)
    check_fname = paste0(check_stub, 
        ".Rda")
    if (file.exists(check_fname)){
        # get the smoothed images
        fdf.run$out_fname = fdf.run$img %>% 
                nii.stub(bn = TRUE) %>% 
                paste0("Reseg_", ., 
                    "_", imod, 
                    "_probability", 
                    "_smoothed") %>% 
                file.path(fdf.run$preddir, .)
        fdf.run$out_fname = gsub("^/dexter",
            "/legacy/dexter",
            fdf.run$out_fname)

        all.p = NULL 
        iimg = 1
        ################################
        # Read smoothed probability image
        # concatenate
        ################################        
        for (iimg in seq(L)){
            out_fname = 
                fdf.run$out_fname[iimg]

            img = readnii(out_fname)
            p = img[ l.keep.ind[[iimg]] ]

            all.p = c(all.p, p)
            print(iimg)
        }
        stopifnot(length(all.p) == 
            nrow(all.df))
        ################################
        # Take the subsampled data
        # So comparable to unsmoothed
        ################################
        test.p = all.p[ !samps ]
        test.y = all.df$Y[ !samps ]
        candidate = mult.df$candidate[!samps]

        ################################
        # Take the candidate data
        # So comparable to unsmoothed
        ################################
        test.p = test.p[candidate]
        test.y = test.y[candidate]

        ################################
        # Performance Metrics
        ################################
        modlist = mod_func(test.p, 
            test.y, 
            fpr.stop = 0.01)
        modlist$mod.perf = NULL

        outname = paste0(check_stub,
            "_smoothed.Rda")
        save(fdf.run, 
            modlist,  
            file = outname)
    }
    print(imod)
}
