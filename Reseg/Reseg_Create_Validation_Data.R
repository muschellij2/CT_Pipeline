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

spec = matrix(c(
    'group', 'c', 1, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

group = opt$group
print(opt)

if (is.null(group)){
    stop("Group is NULL")
}
if (!group %in% c("Test", "Validation")){
    stop("Group is not in Test, Validation")
}

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
    ffdf = fdf
    valid.ind = which(
        ffdf$group == "Validation")
    test.ind = which(
        ffdf$group == "Test")

    # group = "Test"
    if (group == "Validation"){
        subset.ind = valid.ind
    }
    if (group == "Test"){
        subset.ind = test.ind
    }
    fdf.run = ffdf[subset.ind, ]


    moddname = nii.stub(fdf.run$img, 
        bn = TRUE)
    moddname = file.path(fdf.run$outdir, 
        paste0("Reseg_", moddname, 
            "_predictors", 
            adder, ".Rda"))

    all.df = NULL
    L = nrow(fdf.run)
    l.keep.ind = vector(mode="list", 
        length = L)
    imod = 1

   fname = file.path(outdir, 
    paste0("Reseg_Aggregate_data_cutoffs", 
        adder, ".Rda")) 
   load(fname)

    for (imod in seq(L)){
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

        sum(df$Y[!include])/ sum(df$Y)
        sum(include)/ length(include)
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


    sum(df$Y[! df$zval2]) / sum(df$Y)
    sum(1-df$Y[! df$zval2]) / sum(1-df$Y)

    sum(df$Y[! df$zval2_medztemp]) / 
        sum(df$Y)
    sum(1-df$Y[! df$zval2_medztemp]) / 
        sum(1-df$Y)

        df$multiplier = df[, multiplier_col]
        # zval2

        df$candidate = 
            df$multiplier | df$Y == 1

        df = df[ df$candidate, ]

        set.seed(20150215)
        size = min(nrow(df), 1e5)
        samp = sample(nrow(df), size=size)
        df = df[samp, ]

        # l.keep.ind[[imod]] = keep.ind

        df$img = fdf.run$img[imod]
        all.df = rbind(all.df, df)
        rm(list=c("img.pred", "df"))
        print(imod)
    }

    # samp.ind = sample(nrow(df), size= 1e4)
    mult.df = all.df[, 
        !colnames(all.df) %in% keep.colnames]
    all.df = all.df[, 
        colnames(all.df) %in% keep.colnames]


    img = all.df$img
    all.df$img = NULL

    # train = all.df[samps,]
    # test = all.df[!samps,]

    fname = file.path(outdir, 
        paste0("Reseg_", group, 
            "_data", adder, ".Rda"))

    save(all.df, mult.df, img, 
        fdf.run, 
        file = fname)

# }
