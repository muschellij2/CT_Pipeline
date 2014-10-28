####################################################################
## This code is for prediciton of Image Segmentation of CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
####################################################################
####################################################################
rm(list=ls())
library(plyr)

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

#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )

outfile = file.path(outdir, "111_Filenames_with_volumes.Rda")
load(file = outfile)

scen = expand.grid(ttype = c("Affine", "Rigid", "SyN"),
    interpolator = c("Linear", "LanczosWindowedSinc"),
    stringsAsFactors = FALSE)
iscen = 1
addons = paste0(scen$ttype, "_", scen$interpolator) 
addons.diff = paste0(addons, ".diff")
addons.adiff = paste0(addons, ".adiff")

fdf[, addons.adiff] = fdf[, addons.diff] = fdf[, addons] = NA

for (iscen in seq(nrow(scen))){
    
    addon = addons[iscen]
    addon.diff = addons.diff[iscen]
    addon.adiff = addons.adiff[iscen]
    
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


    # template.file = file.path(tempdir, "scct_unsmooth.nii.gz")
    # ss.tempfile = file.path(tempdir, "Skull_Stripped",
    #     "scct_unsmooth_SS_First_Pass_0.1.nii.gz")


    makedir = sapply( fdf$outdir, function(x) {
    	if (!file.exists(x)){
    		dir.create(x, showWarnings =FALSE)
    	}
    })
    fdf$ss = gsub("_Mask", "", fdf$mask)
    iimg = 21

    for (iimg in seq(nrow(fdf))){
        vol = fdf$truevol[iimg]

        x = fdf[iimg,]

        stubfile = function(x, d = NULL, ext = ""){
          b = nii.stub(x, bn=TRUE)
          b = paste0(b, ext)
          file.path(d, b)
        }

        outfile = file.path(x$iddir, "Predictors", 
            paste0("Volume_cutoff_", addon, ".Rda"))

        ####################################
        ## Run both with the Skull Stripped and not skull stripped
        ####################################
        xcut = load(file = outfile) 

        fdf[iimg, addon] = best.cutoff
        fdf[iimg, addon.diff] = vdiff[best.est]
        fdf[iimg, addon.adiff] = adiff[best.est]
        # fdf[iimg, addon.diff] = vdiff[best.est]
        print(iimg)
    }
    print(iscen)
}


cmed = colMedians(as.matrix(fdf[, addons]))
cm = colMeans(as.matrix(fdf[, addons]))
colMins(as.matrix(fdf[, addons]))

(cmed.diff = colMedians(as.matrix(fdf[, addons.diff])))
(cm.diff = colMeans(as.matrix(fdf[, addons.diff])))
(cmax.diff = colMaxs(as.matrix(fdf[, addons.diff])))

(cmed.adiff = colMedians(as.matrix(fdf[, addons.adiff])))
(cm.adiff = colMeans(as.matrix(fdf[, addons.adiff])))
(cmax.adiff = colMaxs(as.matrix(fdf[, addons.adiff])))
