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
roidir = file.path(rootdir, "ROI_data")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

outdir = file.path(basedir, "results")

resdir = file.path(basedir, "results")
load(file=file.path(resdir, "ROI_Filenames.Rda"))

rerun = FALSE
# regs = c("Affine", "Rigid", "SyN")
regs = c("Rigid", "SyN")
# regs = "Rigid"
scen = expand.grid(ttype = regs,
    interpolator = c("Linear"),
    stringsAsFactors = FALSE)
iscen = 2
addons = paste0(scen$ttype, "_", scen$interpolator) 

ids = list.dirs(basedir, recursive=FALSE, full.names=FALSE)
ids = basename(ids)
ids = grep("\\d\\d\\d-(\\d|)\\d\\d\\d", ids, value=TRUE)
length(ids)


iid <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if (is.na(iid)) iid <- 1

all.df$mask = file.path(all.df$iddir, "Skull_Stripped", 
    paste0(nii.stub(all.df$img, bn=TRUE), "_SS_0.01_Mask"))
id = ids[iid]
fdf = all.df[ all.df$id == id, ]

if (nrow(fdf) > 0){


    template.file = file.path(tempdir, "scct_unsmooth.nii.gz")
    ss.tempfile = file.path(tempdir, "Skull_Stripped",
        "scct_unsmooth_SS_0.01.nii.gz")
    ss = readNIfTI(ss.tempfile, reorient=FALSE)


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
        fdf$outdir = file.path(fdf$iddir, outputdir)

        img_ext = paste0("_", ttype, int_ext)

        # CT_Skull_Strip(img = template.file, outfile = ss.tempfile, 
        #     uthresh = 70)
        # ss.tempfile = file.path(tempdir, "Skull_Stripped",
        #     "scct_unsmooth_SS_First_Pass_0.1.nii.gz")

        makedir = sapply( fdf$outdir, function(x) {
        	if (!file.exists(x)){
        		dir.create(x, showWarnings =FALSE)
        	}
        })
        fdf$ss = gsub("_Mask", "", fdf$mask)


        stubfile = function(x, d = NULL, ext = ""){
          b = nii.stub(x, bn=TRUE)
          b = paste0(b, ext)
          file.path(d, b)
        }

        ####################################
        ## Run both with the Skull Stripped and not skull stripped
        ####################################
        iimg = 1

        for (iimg in seq(nrow(fdf))){
            x = fdf[iimg,]


            ofile = stubfile(x$ss, d = x$outdir, ext = img_ext)
            outprefix = ofile

            # outprefix = stubfile(x$ss, d = x$outdir)
            roi.ofile = paste0(ofile, "_ROI.nii.gz")
            mask.ofile = paste0(ofile, "_Mask.nii.gz")
            ofile = paste0(ofile, ".nii.gz")
            print(ofile)
            binary = c(roi.ofile, mask.ofile)
            files = c(ofile, binary)
            ex = all(file.exists(files))
            if (!ex | rerun){
                ss.res = ants_regwrite(
                    filename = x$ss, 
                    correct = FALSE, 
                    other.files = c(x$roi, x$mask), 
                    other.outfiles = c(roi.ofile, mask.ofile),
                    outfile = ofile, retimg = TRUE, 
                    typeofTransform = ttype,
                    template.file = ss.tempfile, 
                    interpolator = interpolator,
                    remove.warp = TRUE,
                    outprefix=NULL)
            }

            if (file.exists(mask.ofile)){
                fslmaths(file= mask.ofile, 
                    outfile = mask.ofile,
                    opts = "-thr 0.5 -bin", 
                    retimg = FALSE)
            }
            if (file.exists(roi.ofile)){
                fslmaths(file= roi.ofile, 
                    outfile = roi.ofile,
                    opts = "-thr 0.5 -bin", 
                    retimg = FALSE) 
            }  
            print(iimg)
        }

    }

}