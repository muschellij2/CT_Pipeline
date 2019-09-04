############################################
## This code is for Image Segmentation of CT
## The code is R but calls out FSL
##
## Author: John Muschelli
## Last updated: May 20, 2014
############################################
############################################
rm(list=ls())
library(plyr)
library(cttools)
library(fslr)
library(mgcv)
library(methods)
homedir = "/Applications"
rootdir = file.path("/Volumes/DATA_LOCAL",
    "Image_Processing")
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = file.path("/dexter/disk2/smart", 
    "stroke_ct", "ident")
}
progdir = file.path(rootdir, "programs")
segdir = file.path(progdir, "Segmentation")

basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")
outdir = file.path(basedir, "results")
rerun = TRUE
#### load voxel data

outfile = file.path(outdir, 
    "Reseg_111_Filenames.Rda")
xxx = load(file = outfile)


fdf$median = fdf$mode = fdf$mean = NA

fdf$thickvol = fdf$zvol = fdf$varslice = 
fdf$gantry = fdf$truevol = NA

# fdf = fdf[c(2,12,17,34,37, 
    # 46, 48, 68, 70,81, 
    #85, 86, 87, 99),]
# dcmtables[, '0018-1152-Exposure']
mods = c("logistic", "lasso", 
    "rf", 
    "gam"
    # , 
    # "cforest"
    )
imod = mods[1]


iimg <- suppressWarnings({
    as.numeric(Sys.getenv("SGE_TASK_ID"))
    })
if (is.na(iimg)) iimg = 34
## 2,12,17,34,37, 46, 48, 68, 
## 70,81, 85, 86, 87, 99
#  has variable slice thickness
## 15 is not correct for sthickness
## 17 & 87 worst - has overlapping 
## slice somewhat
## 75 has all negative
## 71 has no position data
## 13,71,101 has spacing
# for (iimg in seq(nrow(fdf))){
    
    runx = x = fdf[iimg,]
    sortdir = file.path(x$iddir, "Sorted")

    # run_model = function(x, fpr.stop = .1){
    fname = xfname = nii.stub(x$img, bn=TRUE)
    rda = file.path(sortdir, paste0(fname, 
        "_Header_Info.Rda"))
    xrda = load(rda)
    # print(grep("pac", colnames(dcmtables), 
        # value=TRUE))
    # print(iimg)
    for (imod in mods){
        in_stub = x$img %>% 
        nii.stub(bn = TRUE) %>% 
        paste0("Reseg_", ., 
            "_", imod) %>% 
        file.path(x$preddir, .)  
        in_fname = paste0(in_stub, 
            c(
                "_prediction", 
                "_prediction_cc",
             "_smoothed_prediction",
             "_smoothed_prediction_cc"))
        out_fname = paste0(in_fname, 
            "_native.nii.gz")
        out_exists = file.exists(
            out_fname
            )
        
        fname = file.path(x$outdir, 
            paste0("Reseg_", xfname, 
                "_estvolume_", imod, 
                ".Rda"))

        if (
            (all(out_exists)
            & !file.exists(fname)
            ) | rerun
            ){

            n = rep("pred", 
                length = length(out_fname))
            n[grepl("cc", out_fname)] = "cc"
            n[grepl("smooth", out_fname)] = 
                paste0("s", 
                    n[grepl("smooth", 
                        out_fname)]
                    )
            names(out_fname) = n    
            out = llply(out_fname, function(x){
                datatyper(readnii(x) > 0.5)
                }, .progress = "text")

            vols = llply(out, function(x){
                res = get_roi_vol(x, 
                    dcmtables = dcmtables)
                return(res)
                }, .progress = "text")
            est_vols = sapply(vols, `[[`, 
                "truevol")

            save(vols, est_vols,
                file=fname)
        }
    }
    
    print(iimg)
    # print(warnings())
# }
