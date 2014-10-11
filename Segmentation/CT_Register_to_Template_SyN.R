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
library(extrantsr)
library(ANTsR)
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

scen = expand.grid(ttype = c("Affine", "Rigid", "SyN"),
    interpolator = c("Linear", "LanczosWindowedSinc"),
    stringsAsFactors = FALSE)
iscen = 1


# ttype = "Rigid"
ttype = scen$ttype[iscen]
# ttype = "SyN"
# interpolator = c("Linear", "LanczosWindowedSinc")
interpolator = scen$interpolator[iscen]
int_ext = switch(interpolator,
    "Linear" = "",
    "LanczosWindowedSinc" = "_sinc")
outputdir = paste0(ttype, "_Registered")
img_ext = paste0("_", ttype, int_ext)


template.file = file.path(tempdir, "scct_unsmooth.nii.gz")
ss.tempfile = file.path(tempdir, "Skull_Stripped",
    "scct_unsmooth_SS_First_Pass_0.1.nii.gz")
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


makedir = sapply( fdf$outdir, function(x) {
	if (!file.exists(x)){
		dir.create(x, showWarnings =FALSE)
	}
})
fdf$ss = gsub("_Mask", "", fdf$mask)
iimg <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(iimg)) iimg = 108




stubfile = function(x, d = NULL, ext = ""){
  b = nii.stub(x, bn=TRUE)
  b = paste0(b, ext)
  file.path(d, b)
}

####################################
## Run both with the Skull Stripped and not skull stripped
####################################

x = fdf[iimg,]

x$outdir = file.path(x$iddir, outputdir)
if (!file.exists(x$outdir)){
    dir.create(x$outdir, showWarnings =FALSE)
}

ofile = stubfile(x$ss, d = x$outdir, ext = img_ext)
roi.ofile = paste0(ofile, "_ROI.nii.gz")
mask.ofile = paste0(ofile, "_Mask.nii.gz")
ofile = paste0(ofile, ".nii.gz")
print(ofile)
binary = c(roi.ofile, mask.ofile)
files = c(ofile, binary)
ex = all(file.exists(files))
if (!ex){

    ss.res = t1_syn(filename = x$ss, n3correct = FALSE, 
        other.files = c(x$roi, x$mask), 
        other.outfiles = c(roi.ofile, mask.ofile),
        outfile = ofile, retimg = TRUE, 
        typeofTransform = ttype,
        template.file = ss.tempfile, 
        interpolator = interpolator,
        remove.warp = TRUE)
}
if (all(file.exists(binary))){
    fslmaths(file= mask.ofile, 
        outfile = mask.ofile,
        opts = "-thr 0.5 -bin", 
        retimg = FALSE)
    fslmaths(file= roi.ofile, 
        outfile = roi.ofile,
        opts = "-thr 0.5 -bin", 
        retimg = FALSE) 
}  


ofile = stubfile(x$img, d = x$outdir, ext = img_ext)
roi.ofile = paste0(ofile, "_ROI.nii.gz")
mask.ofile = paste0(ofile, "_Mask.nii.gz")
ofile = paste0(ofile, ".nii.gz")
print(ofile)
binary = c(roi.ofile, mask.ofile)
files = c(ofile, binary)
ex = all(file.exists(files))
if (!ex){

    res = t1_syn(filename = x$img, n3correct = FALSE, 
        other.files = c(x$roi, x$mask), 
        other.outfiles = c(roi.ofile, mask.ofile),
        outfile = ofile, retimg = TRUE,
        typeofTransform = ttype,     
        template.file = template.file, 
        interpolator = interpolator,        
        remove.warp = TRUE)
}

if (all(file.exists(binary))){
    fslmaths(file= mask.ofile, 
        outfile = mask.ofile,
        opts = "-thr 0.5 -bin", 
        retimg = FALSE)
    fslmaths(file= roi.ofile, 
        outfile = roi.ofile,
        opts = "-thr 0.5 -bin", 
        retimg = FALSE) 
} 
   