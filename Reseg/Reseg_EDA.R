####################################
## This code is for aggregate models
## CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
####################################
#####################################
rm(list=ls())
library(ggplot2)
library(fslr)
library(mgcv)
library(ROCR)
library(cttools)
set.seed(20150518)
homedir = "/Applications"
rootdir = file.path("/Volumes/DATA_LOCAL", 
    "Image_Processing")
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = file.path("/legacy/dexter/disk2/smart", 
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


    train = train[ order(train$Y), ]
    mydens = function(x){
        g = ggplot(aes(colour = factor(Y)), 
            data = train) + 
        geom_line(stat = "density")
        g = g + aes_string(x = x)
        g
    }

    mylog = function(x, y){
        g = ggplot(aes(y = Y), 
            data = train)
        g = g + aes_string(x = x)
        g = g +
        geom_point(
            position = position_jitter(
                height=0.1),
            size = 0.5
            )
        g = g + geom_smooth(se = FALSE, 
            family = binomial())
        g
    } 

    myscatter = function(x, y){
        g = ggplot(aes(colour = factor(Y)), 
            data = train) + 
        geom_point()
        g = g + aes_string(x = x, y = y)
        g
    }

    mydens("zscore_template")
    mydens("prob_img")
    mydens("moment1")
    mydens("flipped_value")
    mydens("zscore3")

    myscatter("zscore_template", 
        "flipped_value")
    myscatter("dist_centroid", 
        "flipped_value")    
    myscatter("zscore_template", "moment1") 
    myscatter("flipped_value", "moment1") 

    myscatter("moment1", "dist_centroid")


    gam.mod = gam(Y ~ 
        s(zscore_template, flipped_value)
        + s(moment1) 
        + s(dist_centroid)
        , data=train, family= binomial(),
        control = list(trace = TRUE)
        )

    p = predict_gam(gam.mod, newdata = test)

    test.gam.pred = p
    test.gam.pred = as.numeric(test.gam.pred)
    test.gam.pred[ test.gam.pred > 1] = 1
    test$multiplier = 
        test.mult.df$multiplier
    test.gam.pred = test.gam.pred * 
        test$multiplier

    gam.modlist = mod_func(
        test.gam.pred, 
        test$Y, 
        fpr.stop = 0.01)
    gam.modlist$mod = gam.mod
    gam.modlist$mod.time = gam.time
    gam.modlist$mod.perf = NULL
