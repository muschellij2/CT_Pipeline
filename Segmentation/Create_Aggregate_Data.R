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
library(matrixStats)
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

correct = "none"
# options = c("none", "N3", "N4", "N3_SS", "N4_SS",
#         "SyN", "SyN_sinc", "Rigid", "Affine", "Rigid_sinc", 
#         "Affine_sinc")
options = c("none", "N3_SS", "N4_SS", 
      "Rigid", "Rigid_sinc")

#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )

outfile = file.path(outdir, "111_Filenames.Rda")
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
    stopifnot(all(file.exists(outfiles)))
    # fdf = fdf[, ]

    # load(file = file.path(outdir, "Segmentation_Models.Rda"))
    ##############################
    # Run lmod number of models - not all the models - leave out
    ##############################
    fdf.run = fdf[ fdf$group == "Train", ]

    moddname = nii.stub(basename(fdf.run$img))
    moddname = file.path(fdf.run$outdir, 
        paste0(moddname, "_predictors", adder, ".Rda"))

    all.df = NULL
    L = nrow(fdf.run)
    l.keep.ind = vector(mode="list", length = L)
    imod = 1

    for (imod in seq(L)){
        load(moddname[imod])

        df = img.pred$df
        keep.ind = img.pred$keep.ind
        nim = img.pred$nim
        df = df[ keep.ind, ]
        med.ztemp = median(df$zscore_template)

        #df$zval2_medz3 = df$zval2 & (df$zscore3 > med.z3) 
        df$gr_medztemp = (df$zscore_template > med.ztemp) 

        l.keep.ind[[imod]] = keep.ind

        df$img = fdf.run$img[imod]
        all.df = rbind(all.df, df)
        rm(list=c("img.pred", "df"))
        print(imod)
    }

    keep.colnames = colnames(all.df)
    keep.colnames = keep.colnames[ 
        !keep.colnames %in% "gr_medztemp"
        ]

    stopifnot(all(all.df$Y %in% c(0, 1)))
    names(l.keep.ind) = fdf.run$img

    runnames = colnames(all.df)
    nosmooth = c("any_zero_neighbor",
            "thresh", "pct_zero_neighbor")
    runnames = runnames[ !(runnames %in% 
        c("mask", "Y", "img", nosmooth))]


    quants = dlply(all.df[, c(runnames, "Y")], .(Y), 
        function(x) {
        r=t(colQuantiles(as.matrix(x[, runnames]), 
            probs = c(0, 0.001, 0.005, 0.01, 
                0.99, 0.995, 0.999, 1)))
        colnames(r) = runnames
        r
    }, .progress = "text")

    est.cutoffs = quants$`1`[ c("0.5%", "99.5%"), ]

    fname = file.path(outdir, 
        paste0("Aggregate_data_cutoffs", adder, ".Rda"))

    save(est.cutoffs, quants, file = fname)
    
    keepnames = colnames(est.cutoffs)
    include = rep(TRUE, length=nrow(all.df))
    for (icut in keepnames){
        qcuts = est.cutoffs[, icut]
        colname = paste0(icut, ".cutoff")
        all.df[, colname] = all.df[, icut] >= qcuts[1] & 
            all.df[, icut] <= qcuts[2]
        include = include & all.df[, colname]
        print(icut)
    }



    sum(all.df$Y[!include])/ sum(all.df$Y)
    sum(include)/ length(include)
    all.df$include.all = include

    all.df$include = all.df$value >= 30 & all.df$value <= 100


    all.df$zval = all.df[, "zscore3.cutoff"] & all.df$include &
        all.df$pct_thresh.cutoff
    all.df$zval2 = all.df[, "zscore2.cutoff"] & all.df$zval
    all.df$zval_all = all.df[, "zscore_template.cutoff"] & 
        all.df$zval2
    all.df$zval2_medztemp = all.df$zval2 & all.df$gr_medztemp 


    sum(all.df$Y[! all.df$zval2]) / sum(all.df$Y)
    sum(1-all.df$Y[! all.df$zval2]) / sum(1-all.df$Y)

    all.df$multiplier = all.df$zval2_medztemp

    all.df$candidate = all.df$multiplier | all.df$Y == 1

    seed = 20141022
    set.seed(seed)
    ich = which(all.df$Y == 1 & all.df$multiplier) 
    noich = which(all.df$Y != 1 & all.df$multiplier)

    size = 1e5
    prop = .25
    n.ich = ceiling(size*prop)
    n.noich = size - n.ich
    ich.ind = sample(ich, size=n.ich)
    noich.ind = sample(noich, size=n.noich)
    samp.ind = sort(c(ich.ind, noich.ind))

    # samp.ind = sample(nrow(df), size= 1e4)
    samps = seq(nrow(all.df)) %in% samp.ind

    mult.df = all.df[, !colnames(all.df) %in% keep.colnames]
    all.df = all.df[, colnames(all.df) %in% keep.colnames]    


    img = all.df$img
    all.df$img = NULL

    # train = all.df[samps,]
    # test = all.df[!samps,]

    fname = file.path(outdir, 
        paste0("Aggregate_data", adder, ".Rda"))

    save(all.df, mult.df, samps, img, 
        size, prop, n.ich, n.noich, 
        fdf.run, seed, l.keep.ind,
        file = fname)

    train.img = img[samps]
    train = all.df[samps,]
    # tr.img = img[samps,]
    train.mult.df = mult.df[samps,]

    test.img = img[!samps]
    test = all.df[!samps,]
    test.mult.df = mult.df[!samps, ]
    
    test = test[ test.mult.df$candidate, ]
    test.img = test.img[test.mult.df$candidate]
    
    test.mult.df = test.mult.df[ test.mult.df$candidate, ]
    
    fname = file.path(outdir, 
        paste0("Candidate_Aggregate_data", adder, ".Rda"))

    save(train.mult.df, train, 
        test, test.mult.df,
        train.img, test.img,
        # tr.img,
        file = fname)    

# }
