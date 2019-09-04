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


iimg <- as.numeric(
    Sys.getenv("SGE_TASK_ID"))
if (is.na(iimg)) iimg = 28

x = fdf[iimg,]

xmat = load(x$outfile)

df = img.pred$df
keep.ind = img.pred$keep.ind
nim = img.pred$nim

df$mask = df$mask > 0.5
###########################################
# Get cutoffs and make multiplier
###########################################
fname = file.path(outdir, 
    paste0("Reseg_Aggregate_data_cutoffs", 
        adder, ".Rda"))

load(file = fname)

df$gr_medztemp = TRUE
keepnames = colnames(est.cutoffs)
include = rep(TRUE, length=nrow(df))
for (icut in keepnames){
    qcuts = est.cutoffs[, icut]
    colname = paste0(icut, ".cutoff")
    df[, colname] = 
        df[, icut] >= qcuts[1] & 
        df[, icut] <= qcuts[2]
    include = include & df[, colname]
    print(icut)
}
df$include.all = include
df$include = df$value >= 30 & 
df$value <= 100
df$zval = df[, "zscore3.cutoff"] & 
df$include &
    df$pct_thresh.cutoff
df$zval2 = df[, "zscore2.cutoff"] & 
    df$zval
df$zval_all = 
    df[, "zscore_template.cutoff"] & 
    df$zval2
df$zval2_medztemp = df$zval2 & 
    df$gr_medztemp

df$multiplier = df[, multiplier_col]

# mods = c("logistic", "lasso", 
#     "rf", 
#     "gam",
#     "rf_reduced")
mods = "rf_reduced"
# , "gam", "cforest")
imod = mods[1]

for (imod in mods){
    stub_fname = file.path(outdir, 
        paste0("Reseg_Aggregate_models", 
            adder))
    fname = paste0(stub_fname, 
        "_", imod, ".Rda")
    if (file.exists(fname)) {

        load(fname)

        mod = modlist$mod
        cutoff = modlist$mod.dice.coef[1,
            'cutoff']

        f = paste0("predict_", imod)
        func = get(f)
        p = func(mod, df[ df$multiplier, ])
        pimg = remake_img(p, nim, 
            df$multiplier)
        out_fname = x$img %>% 
            nii.stub(bn = TRUE) %>% 
            paste0("Reseg_", ., 
                "_", imod, 
                "_probability") %>% 
            file.path(x$preddir, .)
        writenii(pimg, filename = out_fname)

        mask = remake_img(df$mask, nim)
        pimg = mask_img(pimg, mask)
        # pred = pimg > cutoff
        # cc = spm_bwlabel(pred, k = 100)

        sm.pimg  = mean_image(pimg, 
            nvoxels = 1)
        sm.pimg[abs(sm.pimg) < 
            .Machine$double.eps^0.5 ] = 0
        sm.pimg = niftiarr(nim, sm.pimg)
        sm.pimg[is.na(sm.pimg)]= 0

        sout_fname = paste0(out_fname, 
            "_smoothed")
        writenii(sm.pimg, 
            filename = sout_fname)

        # spred = sm.pimg > cutoff
        # scc = spm_bwlabel(spred, k = 100)

        # les = remake_img(df$Y, nim)
    }
    print(imod)
}


# }