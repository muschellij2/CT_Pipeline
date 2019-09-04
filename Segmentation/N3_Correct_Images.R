#####################################################################
## This code is for N3 Correction of Certain images in CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
#####################################################################
#####################################################################
rm(list=ls())
library(extrantsr)
library(fslr)
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


x = sapply( fdf$outdir, function(x) {
	if (!file.exists(x)){
		dir.create(x, showWarnings =FALSE)
	}
})
x = sapply( fdf$n3dir, function(x) {
    if (!file.exists(x)){
        dir.create(x, showWarnings =FALSE)
    }
})

# irow = 1
irow <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(irow)) irow = 111

# for (irow in seq(nrow(fdf))){
    x = fdf[irow,]

    ss_bias_file = x$n3ssbias_field
    bias_file = x$n3bias_field

    if (!all(file.exists(ss_bias_file, x$n3ssimg))){
        n3ss = bias_correct(file=x$ssimg, outfile = x$n3ssimg, 
            correction = "N3",
            retimg = FALSE,
            reorient = FALSE,
            shrinkfactor = "4",
            "none", "50",  "4", ss_bias_file)
    }
    if (!all(file.exists(bias_file, x$n3img))){
        n3 = bias_correct(file = x$img, outfile = x$n3img, 
            correction = "N3",
            retimg = FALSE,
            reorient = FALSE,
            shrinkfactor = "4",
            "none", "50",  "4", bias_file)
    }
    print(irow)
# }

# for (irow in rev(seq(nrow(fdf)))){
    x = fdf[irow,]
    ss_bias_file = x$n4ssbias_field
    bias_file = x$n4bias_field

    #######################################
# Run N4 on SS image
    #######################################
    if (!all(file.exists(ss_bias_file, x$n4ssimg))){
        args = list(d = 3, 
            i = x$ssimg, 
            o = list(x$n4ssimg, ss_bias_file),
            s = "4")
        N4BiasCorrect_WithField(args)

        args = list(d = 3, 
            i = x$ssimg, 
            o = x$n4ssimg,
            s = "4")
        N4BiasCorrect_WithField(args)
    }

    #######################################
# Run on original image
    #######################################
    # if (!all(file.exists(bias_file, x$n4img))){
        args = list(d = 3, 
            i = x$img, 
            o = list(x$n4img, bias_file),
            s = "4")
        N4BiasCorrect_WithField(args)

        args = list(d = 3, 
            i = x$img, 
            o = x$n4img,
            s = "4")
        N4BiasCorrect_WithField(args)
    # }
    print(irow)
# }

