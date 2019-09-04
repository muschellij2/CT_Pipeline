###########################################################
## This code is for prediction of Image Segmentation of CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
##########################################################
##########################################################
rm(list=ls())
library(plyr)
library(cttools)
library(fslr)
library(ROCR)
library(matrixStats)
library(randomForest)
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

segdir = file.path(progdir, "Segmentation")
source(file.path(segdir, "performance_functions.R"))

outdir = file.path(basedir, "results")

correct = "Rigid"
# options = c("none", "N3", "N4", "N3_SS", "N4_SS",
#         "SyN", "SyN_sinc", "Rigid", 
# "Affine", "Rigid_sinc", 
#         "Affine_sinc")
options = c("none", "N3_SS", "N4_SS",
    "Rigid",  "Rigid_sinc")
# options = c("none", "Rigid_sinc")
# options = "Rigid"

spec = matrix(c(
    'correct', 'c', 1, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

correct = opt$correct
print(opt)



keep.obj = ls()

# for (correct in options){
    
    all.obj = ls()
    rm.obj = all.obj[!(all.obj %in% 
        c(keep.obj, "keep.obj"))]
    rm(list=rm.obj)
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

    #### load voxel data
    outfile = file.path(outdir, "Voxel_Info.Rda")
    load(file=outfile )

    outfile = file.path(outdir, 
        "111_Filenames_with_volumes_stats.Rda")
    load(file = outfile)


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
    outfiles = paste0(outfiles, "_predictors", 
        adder, ".Rda")
    outfiles = file.path(fdf$outdir, outfiles)
    stopifnot(file.exists(outfiles))

    # load(file = file.path(outdir, 
        # "Segmentation_Models.Rda"))
    ##############################
    # Run lmod number of models - 
    # not all the models - leave out
    ##############################
    mod.filename = file.path(outdir, 
        paste0("Collapsed_Models", adder, ".Rda"))
    load(mod.filename)
    vol.sdatas = vol.datas = vol.data
    reses = sreses = res

    cut.vol.data = cut.vol.sdata = cut.vol.tsdata = vol.data
    pauc.cut.vol.data = pauc.cut.vol.sdata = cut.vol.data
    pauc.cut.vol.tsdata = pauc.cut.vol.sdata
    sens.cut.vol.data = sens.cut.vol.sdata = cut.vol.data
    dice.cut.vol.data = dice.cut.vol.sdata = cut.vol.data

    get.pred <- as.numeric(Sys.getenv("SGE_TASK_ID"))
    if (is.na(get.pred)) get.pred = 38
    x = fdf[get.pred,]
    print(get.pred)

# for (get.pred in runpreds){

    iddir = fdf$iddir[get.pred]
    id.outdir = fdf$outdir[get.pred]

    outname = nii.stub(basename(fdf$img[get.pred]))
    outname = file.path(id.outdir, 
        paste0(outname, "_predictors", adder, 
            "_with_subset.Rda"))

    if (!file.exists(outname) | rename){

        predname = nii.stub(basename(fdf$img[get.pred]))
        predname = file.path(id.outdir, 
            paste0(predname, "_predictors", adder, ".Rda"))
        load(predname)
        df = img.pred$df
        stopifnot(all(df$Y %in% c(0, 1)))

        nim = img.pred$nim
        keep.ind = img.pred$keep.ind
        rm(list="img.pred")
        for (i in 1:3) gc() 
        df$include = df$value >= 30 & df$value <= 100
        
        df$mode = fdf$mode[get.pred]

        fname = file.path(outdir, 
            paste0("Aggregate_data_cutoffs", adder, ".Rda"))

        load(file = fname)
        keepnames = colnames(est.cutoffs)
        keepnames = keepnames[!( keepnames %in% 
        c("dist_centroid", "smooth10", "smooth20", "moment2", 
            "moment4", "moment3"))]    

        med.ztemp = median(df$zscore_template[keep.ind])
        df$gr_medztemp = (df$zscore_template > med.ztemp)

        df$skew[is.nan(df$skew)] = 0
        df$kurtosis[is.nan(df$kurtosis)] = 0

        include = rep(TRUE, length=nrow(df))
        for (icut in keepnames){
            qcuts = est.cutoffs[, icut]
            colname = paste0(icut, ".cutoff")
            df[, colname] = df[, icut] >= qcuts[1] & 
                df[, icut] <= qcuts[2]
            include = include & df[, colname]
        }

        df$include.all = include

        df$zval = df[, "zscore3.cutoff"] & df$include &
            df$pct_thresh.cutoff
        df$zval2 = df[, "zscore2.cutoff"] & df$zval
        # df$dist_centroid <= 75
        df$zval_all = df[, "zscore_template.cutoff"] & df$zval2

        df$zval2_medztemp = df$zval2 &  df$gr_medztemp

        # df$subset = df$zval2
        df$subset = df$zval2

        df$in0100 = df$value >= 0 & df$value <= 100
        # df$in20_85 = df$value >= 20 & df$value <= 85
        df$mask = df$mask > 0
        
        ######################################
        # Get volume of ROI
        ######################################  
        vdim = voxdim(nim)
        vres = prod(vdim) / 1000

        Y =  df$Y


        ######################################
        # Keep all ROI = 1, even if not inmask
        ######################################  
        roi.not.in = which(df$Y == 1)
        roi.not.in = roi.not.in[!(roi.not.in %in% keep.ind)]
        keep.ind = sort(c(keep.ind, roi.not.in))


        #### need this because the length of df has changed
        roi.not.in = which(keep.ind %in% roi.not.in)

        df$subset[ roi.not.in ] = FALSE
        df$zval2[ roi.not.in ] = FALSE
        df$zval[ roi.not.in ] = FALSE
        df$zval_all[ roi.not.in ] = FALSE
        df$zval2_medztemp[ roi.not.in ] = FALSE
        df$include[ roi.not.in ] = FALSE
        df$include.all[ roi.not.in ] = FALSE

        ddf = df[keep.ind,]

        ddf = ddf[, !grepl("[.]cutoff", colnames(ddf))]
        rownames(ddf) = NULL


        save(ddf, 
            nim, 
            keep.ind, 
            file = outname)

    }