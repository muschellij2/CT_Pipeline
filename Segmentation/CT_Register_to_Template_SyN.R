###################################################################
## This code is for prediciton of Image Segmentation of CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
####################################################################
####################################################################
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

#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )

outfile = file.path(outdir, "111_Filenames_with_volumes.Rda")
load(file = outfile)


# regs = c("Affine", "Rigid", "SyN")
regs = "Rigid"
scen = expand.grid(ttype = regs,
    interpolator = c("Linear", "LanczosWindowedSinc"),
    stringsAsFactors = FALSE)
iscen = 1

# scen = scen[2,, drop=FALSE]
addons = paste0(scen$ttype, "_", scen$interpolator) 

fdf[, addons] = NA

for (iscen in seq(nrow(scen))){
    
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


    template.file = file.path(tempdir, "scct_unsmooth.nii.gz")
    ss.tempfile = file.path(tempdir, "Skull_Stripped",
        "scct_unsmooth_SS_First_Pass_0.1.nii.gz")

    makedir = sapply( fdf$outdir, function(x) {
    	if (!file.exists(x)){
    		dir.create(x, showWarnings =FALSE)
    	}
    })
    fdf$ss = gsub("_Mask", "", fdf$mask)
    iimg <- as.numeric(Sys.getenv("SGE_TASK_ID"))
    if (is.na(iimg)) iimg = 21
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

    ofile = stubfile(x$ss, d = x$outdir, ext = img_ext)
    # outprefix = stubfile(x$ss, d = x$outdir)
    zfile = file.path( x$preddir, 
        paste0(nii.stub(x$ss, bn=TRUE), 
            "_template_zscore.nii.gz" ))
    zofile = stubfile(zfile, d = x$outdir, ext = img_ext)

    outprefix = ofile

    roi.ofile = paste0(ofile, "_ROI.nii.gz")
    mask.ofile = paste0(ofile, "_Mask.nii.gz")
    ofile = paste0(ofile, ".nii.gz")
    zofile = paste0(zofile, ".nii.gz")
    print(ofile)
    binary = c(roi.ofile, mask.ofile)
    files = c(ofile, binary)
    ex = all(file.exists(files))
    # if (!ex){
        ss.res = ants_regwrite(filename = x$ss, 

            correct = FALSE, 
            other.files = c(x$roi, x$mask, zfile), 
            other.outfiles = c(roi.ofile, mask.ofile, zofile),
            outfile = ofile, retimg = TRUE, 
            typeofTransform = ttype,
            template.file = ss.tempfile, 
            interpolator = interpolator,
            remove.warp = TRUE,
            outprefix=outprefix)
        
        # roi =readNIfTI(roi.ofile, 
        #     reorient=FALSE)
        # vals = roi[roi > 0]
        # cuts = seq(0.01, 1, by=0.001)
        # vres = prod(pixdim(roi)[2:4])
        # stopifnot(vres == 1)
        # est = sapply(cuts, function(x) sum(vals >= x)/1000)
        # vdiff = est - vol
        # adiff = abs(est - vol)
        # best.est = which.min(adiff)
        # best.cutoff = cuts[best.est]
        # fdf[iimg, addon] = best.cutoff
        # outfile = file.path(x$iddir, "Predictors", 
        #     paste0("Volume_cutoff_", addon, ".Rda"))
        # save(best.cutoff, vdiff, adiff, cuts, 
        #     est, best.est, file = outfile) 
    # }
    # if (file.exists(mask.ofile)){
    #     fslmaths(file= mask.ofile, 
    #         outfile = mask.ofile,
    #         opts = "-thr 0.5 -bin", 
    #         retimg = FALSE)
    # }
    # if (file.exists(roi.ofile)){
    #     fslmaths(file= roi.ofile, 
    #         outfile = roi.ofile,
    #         opts = "-thr 0.5 -bin", 
    #         retimg = FALSE) 
    # }  


    #########################################        
    # Non-skull stripped data
    #########################################
    # ofile = stubfile(x$img, d = x$outdir, ext = img_ext)
    # roi.ofile = paste0(ofile, "_ROI.nii.gz")
    # mask.ofile = paste0(ofile, "_Mask.nii.gz")
    # ofile = paste0(ofile, ".nii.gz")
    # print(ofile)
    # binary = c(roi.ofile, mask.ofile)
    # files = c(ofile, binary)
    # ex = all(file.exists(files))
    # if (!ex){

    #     res = ants_regwrite(filename = x$img, 
    #         n3correct = FALSE, 
    #         other.files = c(x$roi, x$mask), 
    #         other.outfiles = c(roi.ofile, mask.ofile),
    #         outfile = ofile, retimg = TRUE,
    #         typeofTransform = ttype,     
    #         template.file = template.file, 
    #         interpolator = interpolator,        
    #         remove.warp = TRUE)
    #     roi =readNIfTI(roi.ofile, 
    #         reorient=FALSE)
    #     vals = roi[roi > 0]
    #     cuts = seq(0.01, 1, by=0.001)
    #     vres = prod(pixdim(roi)[2:4])
    #     stopifnot(vres == 1)
    #     est = sapply(cuts, function(x) sum(vals > x)/1000)
    #     vdiff = abs(est - vol)
    #     best.est = which.min(vdiff)
    #     best.cutoff = cuts[best.est]
    #     fdf$ssbest.cut[iimg] = best.cutoff
    #
    # }

    # if (file.exists(mask.ofile)){
    #     fslmaths(file= mask.ofile, 
    #         outfile = mask.ofile,
    #         opts = "-thr 0.5 -bin", 
    #         retimg = FALSE)
    # }
    # if (file.exists(roi.ofile)){    
    #     fslmaths(file= roi.ofile, 
    #         outfile = roi.ofile,
    #         opts = "-thr 0.5 -bin", 
    #         retimg = FALSE) 
    # } 
    print(iscen)
}

