#####################################################################
## This code is for prediciton of Image Segmentation of CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
#####################################################################
#####################################################################
rm(list=ls())
library(plyr)
library(cttools)
library(fslr)
library(ROCR)
library(matrixStats)
library(ggplot2)
library(smallPDF)
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

correct = "N3"
options = c("none", "N3", "N4", "N3_SS", "N4_SS",
        "SyN", "SyN_sinc", "Rigid", "Affine")
for (correct in options){
    
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


    short_predict = function(object, newdata, 
    	lthresh=  .Machine$double.eps^0.5){
    	tt <- terms(object)
    	Terms <- delete.response(tt)
        m <- model.frame(Terms, newdata, 
        	na.action = na.pass, 
            xlev = object$xlevels)
        if (is.null(cl <- attr(Terms, "dataClasses"))) 
                stop("no dataclasses")
        X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
        # p <- object$rank
        beta <- object$coefficients
        beta = beta[ !is.na(beta) ]
        predictor = drop(X[, names(beta), drop=FALSE ] %*% beta)
       
    	predictor <- family(object)$linkinv(predictor)
    	predictor[ predictor < lthresh] = 0
    	predictor
    }

    #### load voxel data
    outfile = file.path(outdir, "Voxel_Info.Rda")
    load(file=outfile )

    outfile = file.path(outdir, "111_Filenames.Rda")
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
    outfiles = paste0(outfiles, "_predictors.Rda")
    outfiles = file.path(fdf$outdir, outfiles)
    fdf = fdf[file.exists(outfiles), ]

    # load(file = file.path(outdir, "Segmentation_Models.Rda"))
    ##############################
    # Run lmod number of models - not all the models - leave out
    ##############################
    lmod = 10
    fdf.run = fdf[seq(lmod), ]

    moddname = nii.stub(basename(fdf.run$img))
    moddname = file.path(fdf.run$outdir, 
        paste0(moddname, "_predictors", adder, ".Rda"))

    all.df = NULL
    for (imod in seq(lmod)){
        load(moddname[imod])

        df = img.pred$df
        keep.ind = img.pred$keep.ind
        nim = img.pred$nim
        df = df[ keep.ind, ]

        df$img = fdf.run$img[imod]
        all.df = rbind(all.df, df)
        rm(list=c("img.pred", "df"))
        print(imod)
    }


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
    train = all.df[samps,]
    test = all.df[!samps,]

    runnames = names(train)
    nosmooth = c("any_zero_neighbor",
            "thresh", "pct_zero_neighbor")
    runnames = runnames[ !(runnames %in% 
        c("mask", "Y", "img", nosmooth))]
    pdfname = file.path(outdir, 
        paste0("Aggregate_Data_Plots", adder, ".pdf"))
    # pdf(pdfname)
    pdfobj = smallpdf()
        g = ggplot(train, aes(y = Y)) + geom_point(
            position = position_jitter(w = 0.0, h = 0.2))
        g2 = g + stat_smooth(se=FALSE, method = "gam", 
            formula = y ~ s(x, bs = "cs"),
            family = binomial())   
        for (iname in runnames){
            gg = g2 + aes_string(x = iname)
            gg = gg + ggtitle(iname)
            print(gg)
            print(iname)
        }
        for (iname in nosmooth){
            gg = g + aes_string(x = iname)
            gg = gg + ggtitle(iname)
            print(gg)
            print(iname)
        }
    smallpdf.off(pdfname = pdfname, 
        pdfobj = pdfobj, 
        clean = TRUE)

}