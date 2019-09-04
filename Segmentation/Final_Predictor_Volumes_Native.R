##################################################################
## This code is for prediciton of Image Segmentation of CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
##################################################################
##################################################################
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

outfile = file.path(outdir, 
    "111_Filenames_with_volumes_stats.Rda")
load(file = outfile)


regs = "Rigid"
scen = expand.grid(ttype = regs,
    interpolator = c("Linear"),
    stringsAsFactors = FALSE)
iscen = 1
# scen = scen[2,, drop=FALSE]
addons = paste0(scen$ttype, "_", scen$interpolator) 

fdf[, addons] = NA

types = c("_zval2", '_zval2_medztemp')
# , "_zval2"
# "_include_all", 
type = types[1]

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

fdf$hdr = file.path(fdf$iddir, "Sorted",
    paste0(nii.stub(fdf$img, bn=TRUE), 
        "_Header_Info.Rda"))

iimg <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(iimg)) iimg = 3

fdf$vol.raw = fdf$vol.pred = fdf$vol.pred_0.5 = NA


for (iimg in seq(nrow(fdf))){

    x = fdf[iimg,]    
    vol = x$truevol

    stubfile = function(x, d = NULL, ext = ""){
      b = nii.stub(x, bn=TRUE)
      b = paste0(b, ext)
      file.path(d, b)
    }

    ####################################
    ## Run both with the Skull Stripped and not skull stripped
    ####################################
    x$preddir = x$outdir
    x$outdir = file.path(x$iddir, outputdir)
    hdrload = load(x$hdr)



    # cn = c("mod_agg", "gam")
    cn = "mod_agg"
    outimg = nii.stub(x$img, bn=TRUE)
    outimg = file.path(x$preddir, 
        paste0(outimg, "_", cn, adder, type, "_prediction"))
    outfile = paste0(outimg, "_native.nii.gz")

    iout = 1
    # for (iout in seq_along(outfile)){

        native.roi = cal_img(readNIfTI(outfile[iout], 
            reorient=FALSE))
        native.pred = cal_img(native.roi > 0)
        native.pred_0.5 = cal_img(native.roi > 0.5)


        vol.raw = get_roi_vol(native.roi, dcmtables)$truevol
        vol.pred = get_roi_vol(native.pred, dcmtables)$truevol
        vol.pred_0.5 = get_roi_vol(native.pred_0.5, 
            dcmtables)$truevol
        print(outfile[iout])
        print(vol)
        print(vol.raw)
        print(vol.pred)
        print(vol.pred_0.5)
        fdf$vol.raw[iimg] = vol.raw
        fdf$vol.pred[iimg] = vol.pred
        fdf$vol.pred_0.5[iimg] = vol.pred_0.5
    # }
    # 
    print(iimg)
}
    
    # print(iscen)
# }

outfile = file.path(outdir, 
    paste0("111_Filenames_with_volumes", "_stats", 
        "_predictions", ".Rda"))
save(fdf, file = outfile)