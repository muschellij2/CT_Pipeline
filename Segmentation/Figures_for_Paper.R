###########################################################
## This calculates the volumes for the Rigid back trans
## formed data
##
## Author: John Muschelli
## Last updated: May 20, 2014
###########################################################
###########################################################
rm(list=ls())
library(plyr)
library(cttools)
library(fslr)
library(ROCR)
library(matrixStats)
library(extrantsr)
library(ANTsR)
library(scales)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
segdir = file.path(progdir, "Segmentation")

basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

outdir = file.path(basedir, "results")

cn = c("mod_agg")
icut = c("_dice")
types = "_zval2"
type = types[1]
correct = "none"
adder = ""


scut.filename = file.path(outdir, 
paste0("Smooth_Model_Cutoffs", adder, ".Rda"))

load(file=scut.filename)

#########
# Loading in fdf
#########
outfile = file.path(outdir, 
    "111_Filenames_with_volumes_stats.Rda")
load(file = outfile)
group = "Test"


ss_smooth = c(outer( paste0(cn, adder, 
    "_smoothed", adder), 
    type, paste0))

fdf$outimg = nii.stub(fdf$img, bn=TRUE)
fdf$outimg = file.path(fdf$outdir ,
    paste0(fdf$outimg, "_", ss_smooth))
fdf$outfile = paste0(fdf$outimg, ".nii.gz")         

fdf.run = fdf[ fdf$group == group, ]

#####################
# Loading in DICE
#####################
outfile = file.path(outdir, 
    paste0("Model_performance_results", adder, type, 
        ".Rda")
    )
xx = load(outfile)
dice.mod.sdice = dice.mod.sdice[ fdf$group == group, 
    cn, drop = FALSE]
ranks = rank(dice.mod.sdice)
inds = floor(quantile(1:nrow(fdf.run)))
pick = which(ranks %in% inds)
pick = pick[ order(ranks[pick])]

fdf.pick = fdf.run[pick, ]
dices = dice.mod.sdice[pick, , drop=FALSE]
# types = c("_zval2", '_zval2_medztemp')

fdf = fdf.pick
iimg = 4

for (iimg in seq(nrow(fdf))){

    q = names(inds)[iimg]
    q = as.numeric(gsub("%", "", q)) 
    q = sprintf("%03.0f", q)

    pngname = file.path(outdir, 
        paste0("Figure_DSI_Quantile_", q, ".png"))
    if (!file.exists(pngname)){

        x = fdf[iimg,]
        ####################################
        # Read in ROI and cross reference it with Prediction
        ####################################
        x$preddir = x$outdir

        img = readNIfTI(x$ssimg, 
                reorient=FALSE)

        roi = readNIfTI(x$roi, 
                reorient=FALSE)
        
        # ortho2(img, roi, xyz = xyz)

        native.roi = readNIfTI(x$outfile, 
                reorient=FALSE)
        cutval = all.sdice.cuts[cn]
        native.pred = native.roi > cutval

        ###########################
        ### Drop empty image dimensions
        ###########################
        o = dropEmptyImageDimensions(img, 
            other.imgs = list(roi, 
                native.roi, native.pred))
        img = o$outimg
        roi = o$other.imgs[[1]]
        native.roi = o$other.imgs[[2]]
        native.pred = o$other.imgs[[3]]

        # ortho2(img, native.pred, xyz = xyz)

        xyz=xyz(roi)
        xyz[3] = which.max(apply(roi, 3, sum))


        cols = c("#56B4E9", "#D55E00", "#009E73")

        diff = niftiarr(img, NA)
        # false negative
        diff[ roi == 1 & native.pred == 0] = 1
        # false positive
        diff[ roi == 0 & native.pred == 1] = 2
        # true positive
        diff[ roi == 1 & native.pred == 1] = 3
        diff = cal_img(diff)


        plevs = c("False Negative", 
            "False Positive", 
            "True Positive")

        png(pngname,
            res = 600, units = "in", height=7, 
            width = 7, type= "cairo")
        ortho2(img, diff, 
            window=c(0, 100),
            # don't do alpha blending
            col.y = cols,
            xyz = xyz, 
            addlegend = TRUE,
            legend=plevs, 
            leg.col=cols, 
            leg.cex=1.5,
            ybreaks = c(0, 1.1, 2.1, 3.1))
        dev.off()
    }
    print(iimg)
}

