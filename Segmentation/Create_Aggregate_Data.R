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

correct = "SyN"
options = c("none", "N3", "N4", "N3_SS", "N4_SS",
        "SyN", "SyN_sinc", "Rigid", "Affine")

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
        "Affine" = "_Affine")

    filename = file.path(outdir, 
        paste0("Collapsed_Models", adder, ".Rda"))
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
    fdf = fdf[file.exists(outfiles), ]

    # load(file = file.path(outdir, "Segmentation_Models.Rda"))
    ##############################
    # Run lmod number of models - not all the models - leave out
    ##############################
    fdf.run = fdf[run.ind, ]

    moddname = nii.stub(basename(fdf.run$img))
    moddname = file.path(fdf.run$outdir, 
        paste0(moddname, "_predictors", adder, ".Rda"))

    all.df = NULL
    L = nrow(fdf.run)
    l.keep.ind = vector(mode="list", length = L)
    for (imod in seq(L)){
        load(moddname[imod])

        df = img.pred$df
        keep.ind = img.pred$keep.ind
        nim = img.pred$nim
        df = df[ keep.ind, ]
        l.keep.ind[[imod]] = keep.ind

        df$img = fdf.run$img[imod]
        all.df = rbind(all.df, df)
        rm(list=c("img.pred", "df"))
        print(imod)
    }

    names(l.keep.ind) = fdf.run$img

    runnames = colnames(all.df)
    nosmooth = c("any_zero_neighbor",
            "thresh", "pct_zero_neighbor")
    runnames = runnames[ !(runnames %in% 
        c("mask", "Y", "img", nosmooth))]


    quants = dlply(all.df[, c(runnames, "Y")], .(Y), function(x) {
        r=t(colQuantiles(x[, runnames], 
            probs = c(0, 0.001, 0.01, 0.99, 0.999, 1)))
        colnames(r) = runnames
        r
    })

    est.cutoffs = quants$`1`[ c("0.1%", "99.9%"), ]

    fname = file.path(outdir, 
        paste0("Aggregate_data_cutoffs", adder, ".Rda"))

    save(est.cutoffs, quants, file = fname)
    
    seed = 20141022
    set.seed(seed)
    ich = which(all.df$Y == 1)
    noich = which(all.df$Y != 1)

    size = 1e5
    prop = .25
    n.ich = ceiling(size*prop)
    n.noich = size - n.ich
    ich.ind = sample(ich, size=n.ich)
    noich.ind = sample(noich, size=n.noich)
    samp.ind = sort(c(ich.ind, noich.ind))

    # samp.ind = sample(nrow(df), size= 1e4)
    samps = seq(nrow(all.df)) %in% samp.ind
    
    img = all.df$img
    all.df$img = NULL    

    # train = all.df[samps,]
    # test = all.df[!samps,]

    fname = file.path(outdir, 
        paste0("Aggregate_data", adder, ".Rda"))

    save(all.df, samps, img, 
        size, prop, n.ich, n.noich, 
        fdf.run, seed, l.keep.ind,
        file = fname)

# }
