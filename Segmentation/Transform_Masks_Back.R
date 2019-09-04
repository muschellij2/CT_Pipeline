###################################################################
## This code is for prediciton of Image Segmentation of CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
###################################################################
###################################################################
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

outdir = file.path(basedir, "results")

template.file = file.path(tempdir, "scct_unsmooth.nii.gz")
ss.tempfile = file.path(tempdir, "Skull_Stripped",
    "scct_unsmooth_SS_First_Pass_0.1.nii.gz")

#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )

outfile = file.path(outdir, "111_Filenames_with_volumes.Rda")
load(file = outfile)

# types = c("_zval2", '_zval2_medztemp')
types = c("_zval2")
# , "_zval2"
# "_include_all", 
type = types[1]

regs = "Rigid"
scen = expand.grid(ttype = regs,
    interpolator = c("Linear"),
    stringsAsFactors = FALSE)
iscen = 1
# scen = scen[2,, drop=FALSE]
addons = paste0(scen$ttype, "_", scen$interpolator) 

fdf[, addons] = NA

# for (iscen in seq(nrow(scen))){
    
    addon = addons[iscen]
    
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
    correct = gsub("^_", "", img_ext)
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

    makedir = sapply( fdf$outdir, function(x) {
    	if (!file.exists(x)){
    		dir.create(x, showWarnings =FALSE)
    	}
    })
    fdf$ss = gsub("_Mask", "", fdf$mask)
    iimg <- as.numeric(Sys.getenv("SGE_TASK_ID"))
    if (is.na(iimg)) iimg = 38
    vol = fdf$truevol[iimg]

    stubfile = function(x, d = NULL, ext = ""){
      b = nii.stub(x, bn=TRUE)
      b = paste0(b, ext)
      file.path(d, b)
    }

    ####################################
    ## Run both with the Skull Stripped and not skull stripped
    ####################################

    x = fdf[iimg,]
    x$preddir = x$outdir
    x$outdir = file.path(x$iddir, outputdir)
    if (!file.exists(x$outdir)){
        dir.create(x$outdir, showWarnings =FALSE)
    }

    cn = c("mod_agg", "gam", "rf")
    ss = c(outer( paste0(cn, adder), types, paste0))
    ss_smooth = c(outer( paste0(cn, adder, "_smoothed", adder), 
        types, paste0))

    xout = outimg = nii.stub(x$img, bn=TRUE)
    outimg = file.path(x$preddir, 
        paste0(outimg, "_", ss, "_prediction"))
    xout = file.path(x$preddir, 
        paste0(xout, "_", ss_smooth))
    outimg = c(outimg, xout)


    outfile = paste0(outimg, "_native.nii.gz")
    outimg = paste0(outimg, '.nii.gz')
    #     
    # for (iout in seq_along(outfile)){
    iout = seq_along(outimg)
    
        # template = readNIfTI(ss.tempfile, reorient=FALSE)
        # roi = antsImageRead(outimg[iout], 3)

        ss.res = ants_regwrite(filename = ss.tempfile,
            correct = FALSE, 
            other.files = c(outimg[iout]), 
            other.outfiles = c(outfile[iout]),
            outfile = paste0(tempfile(), ".nii.gz"), 
            retimg = FALSE, 
            typeofTransform = ttype,
            template.file = x$ss, 
            interpolator = interpolator,
            remove.warp = TRUE,
            outprefix=NULL)

    # }

    # native.roi = readNIfTI(outfile, reorient=FALSE)
    # native.pred = cal_img(native.roi > 0)
    # ssimg = readNIfTI(x$ss, reorient=FALSE)

    
    # print(iscen)
# }

