#######################################
## This code is for aggregate prediction 
## of Image Segmentation of 
## CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
########################################
rm(list=ls())
library(methods)
library(plyr)
library(cttools)
library(fslr)
library(matrixStats)
library(getopt)
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

correct = "Rigid"
options = c("Rigid")

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
            dir.create(x, 
                showWarnings =FALSE)
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
    stopifnot(all(file.exists(outfiles)))
    # fdf = fdf[file.exists(outfiles), ]

    # load(file = file.path(outdir, 
        # "Segmentation_Models.Rda"))
    ##############################
    # Run lmod number of models - 
    # not all the models - leave out
    ##############################


    moddname = nii.stub(fdf$img, 
        bn = TRUE)
    moddname = file.path(fdf$outdir, 
        paste0("Reseg_", moddname, 
            "_predictors", 
            adder, ".Rda"))
    imod = 1

    fname = file.path(outdir, 
        paste0(
            "Reseg_Aggregate_data_cutoffs", 
            adder, ".Rda")) 
    load(fname)
    fdf$includeY1 = NA
    fdf$includeY0 = NA

    for (imod in seq(nrow(fdf))){
        load(moddname[imod])

        df = img.pred$df
        keep.ind = img.pred$keep.ind
        nim = img.pred$nim
        df$ind = seq(nrow(df))
        df = df[ keep.ind, ]
        med.ztemp = median(df$zscore_template)

        df$gr_medztemp = 
            (df$zscore_template > med.ztemp) 

        keep.colnames = colnames(df)

        stopifnot(all(df$Y %in% c(0, 1)))
  
        keepnames = colnames(est.cutoffs)
        include = rep(TRUE, length=nrow(df))
        for (icut in keepnames){
            qcuts = est.cutoffs[, icut]
            colname = paste0(icut, ".cutoff")
            df[, colname] = 
                df[, icut] >= qcuts[1] & 
                df[, icut] <= qcuts[2]
            include = include & df[, colname]
            # print(icut)
        }

        df$include.all = include

        df$include = 
            df$value >= 30 & df$value <= 100


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

        mult = df[, multiplier_col]

        fdf$includeY1[imod] = 
            sum(df$Y[mult]) / sum(df$Y)
        fdf$includeY0[imod] = 
            sum(1-df$Y[mult]) / sum(1-df$Y)
        print(imod)
    }

#### load voxel data
outrda = file.path(outdir, 
    "Reseg_111_Filenames_with_Exclusions.Rda")
save(fdf, file = outrda)