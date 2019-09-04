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
library(extrantsr)
library(methods)
library(magrittr)
rerun = TRUE
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
outdir = file.path(basedir, 
    "results")

segdir = file.path(progdir, 
    "Reseg")
source(file.path(segdir, 
    "Reseg_performance_functions.R"))

correct = "Rigid"

options = c("Rigid")
interpolator = "Linear"
# options = c("none", "Rigid")

#### load voxel data
outfile = file.path(outdir, 
    "Reseg_111_Filenames.Rda")
load(file = outfile)

for (i in 1:3) gc()
correct = match.arg(correct, options)
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

##############################
# Keeping files where predictors exist
##############################
iimg <- as.numeric(
    Sys.getenv("SGE_TASK_ID"))
if (is.na(iimg)) iimg = 21

x = fdf[iimg,]

stubfile = function(x, d = NULL, ext = ""){
  b = nii.stub(x, bn=TRUE)
  b = paste0(b, ext)
  file.path(d, b)
}
    
int_ext = switch(interpolator,
    "Linear" = "",
    "LanczosWindowedSinc" = "_sinc")
outputdir = paste0(correct, "_Registered")
img_ext = paste0("_", correct, int_ext)

x = fdf[iimg,]
x$preddir = x$outdir
x$outdir = file.path(x$iddir, outputdir)


ofile = stubfile(x$ss, 
    d = x$outdir, 
    ext = img_ext)
outprefix = ofile
omat = paste0(outprefix, 
    "0GenericAffine.mat")

###########################################
# Get cutoffs and make multiplier
###########################################
# mods = c("logistic", "lasso", 
#     "rf", 
#     "gam", 
#     "cforest")
mods = "rf_reduced"
# , "gam", "cforest")
imod = mods[1]

fixed = x$ss
roi = readnii(x$roi)
mask = readnii(x$mask)
mask = mask > 0 | roi > 0
vals = roi[mask]
vdim = prod(voxdim(roi))


for (imod in mods){
    
    in_stub = x$img %>% 
        nii.stub(bn = TRUE) %>% 
        paste0("Reseg_", ., 
            "_", imod) %>% 
        file.path(x$preddir, .)             
    
    in_fname = paste0(in_stub, 
            c(
                "_probability", 
                "_probability_smoothed",
                "_prediction", 
                "_prediction_cc",
             "_smoothed_prediction",
             "_smoothed_prediction_cc"))
    in_fname_full = paste0(in_fname, 
        ".nii.gz")
    check_exists = file.exists(
        in_fname_full
        )
    out_fname = paste0(in_fname, 
        "_native.nii.gz")
    out_exists = file.exists(
        out_fname
        )

    if (all(check_exists)){
        if (!all(out_exists) |
            rerun ){
            message("# transforming data\n")

            out = ants_apply_transforms(
            fixed = fixed, 
            moving = in_fname_full,
            transformlist = omat, 
            whichtoinvert = 1)

            mapply(function(x, fname){
                writenii(x, filename = fname)
            }, out, out_fname)
        } else {
            message("# reading in data\n")
            out = lapply(out_fname, readnii)
        }
    } else {
        print("Images do not exist!")
    }

    keep_preds = grepl("prediction", 
        out_fname)
    out_preds = out[keep_preds]

    results = lapply(out_preds, 
        function(x){
        e = extrantsr:::sim(
            vals,
            x[mask] > 0.5)
        e$truevol = e$truevol * vdim
        e$estvol = e$estvol * vdim
        e
    })
    out_pred_names = out_fname[keep_preds]
    n = rep("pred", 
        length = length(out_pred_names))
    n[grepl("cc", out_pred_names)] = "cc"
    n[grepl("smooth", out_pred_names)] = 
        paste0("s", 
            n[grepl("smooth", out_pred_names)]
            )
    names(results) = n

    res_rda = paste0(x$res_rda_stub, 
        "_", imod, "_native.Rda")
    save(results, file = res_rda)
    print(imod)
}




