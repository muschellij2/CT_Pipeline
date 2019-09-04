############################################
## This calculates the volumes for the 
# Rigid back trans
## formed data
##
## Author: John Muschelli
## Last updated: 2015/11/10
#############################################
#############################################
rm(list=ls())
library(plyr)
library(dplyr)
library(cttools)
library(fslr)
library(ROCR)
library(matrixStats)
library(extrantsr)
library(ANTsR)
library(scales)
library(knitr)
homedir = "/Applications"
rootdir = 
    "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = 
    "/legacy/dexter/disk2/smart/stroke_ct/ident"
}

basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

outdir = file.path(basedir, "results")
rerun = TRUE


#########
# Loading in fdf
#########
outfile = file.path(outdir, 
    paste0("Reseg_Results.Rda"))
load(file = outfile)
run_group = c("Test", "Validation")

fdf.run = fdf[ fdf$group %in% run_group, ]
# need to subset long
if (!"group" %in% colnames(long)){
  long = merge(long, fdf.run[, c("id", "group")], all.x = TRUE)
}
long = long[ long$group %in% run_group, ]


#####################
# Loading in DICE
#####################
run_cutoff = "scc"
run_mod = "rf"
correct = "Native"
iapp = switch(correct,
    Native = "_native", 
    Rigid = "")
dice = filter(long, 
    mod %in% run_mod, 
    app %in% correct,
    group %in% run_group,
    cutoff %in% run_cutoff)
dice = select(dice, iimg, cutoff, 
    dice, group, app)

iimg = 4

qs = quantile(dice$dice)
ranks = rank(dice$dice)
inds = floor(quantile(1:nrow(fdf.run)))
pick = which(ranks %in% inds)
pick = pick[ order(ranks[pick])]

fdf.run$preddir = gsub("^/dexter",
    "/legacy/dexter",
    fdf.run$preddir )
fdf.run$img = gsub("^/dexter",
    "/legacy/dexter",
    fdf.run$img )
fdf.run$rig_ssimg = gsub("^/dexter",
    "/legacy/dexter",
    fdf.run$rig_ssimg )
fdf.run$ssimg = gsub("^/dexter",
    "/legacy/dexter",
    fdf.run$ssimg )

fdf.run$roi = gsub("^/dexter",
    "/legacy/dexter",
    fdf.run$roi )
fdf.run$rig_ssroi = gsub("^/dexter",
    "/legacy/dexter",
    fdf.run$rig_ssroi )

fdf.pick = fdf.run[pick, ]

dice = dice[pick, , drop=FALSE]

iimg = 2

for (iimg in seq(nrow(fdf.pick))){

    q = names(inds)[iimg]
    q = as.numeric(gsub("%", "", q)) 
    q = sprintf("%03.0f", q)

    pngname = file.path(outdir, 
        paste0("Reseg_Figure_DSI_Quantile_", 
            q, iapp, ".png"))

    max_spngname = file.path(outdir, 
        paste0("Reseg_Slice_", "Max", 
            "_DSI_Quantile_", 
            q, iapp, ".png"))    
    med_spngname = file.path(outdir, 
        paste0("Reseg_Slice_", "Med", 
            "_DSI_Quantile_", 
            q, iapp, ".png"))  
    if (!all(
        file.exists(pngname),
        file.exists(max_spngname),
        file.exists(med_spngname)
        ) | rerun
    ) {

        xx = x = fdf.pick[iimg,]
        run_dice = dice$dice[iimg]    

        in_stub = x$img %>% 
                nii.stub(bn = TRUE) %>% 
                paste0("Reseg_", ., 
                    "_", run_mod) %>% 
                file.path(x$preddir, .)
        pred = paste0(in_stub, 
                "_prediction", iapp)
        cc = paste0(in_stub, 
                "_prediction_cc", iapp)
        sin_fname = paste0(in_stub, 
            "_smoothed")
        spred = paste0(sin_fname, 
                "_prediction", iapp)
        scc = paste0(sin_fname, 
                "_prediction_cc", iapp)  

        img_name = switch(run_cutoff,
            cc = cc,
            pred = pred,
            spred = spred,
            scc = scc)
        ####################################
        # Read in ROI and cross reference 
        # it with Prediction
        ####################################
        x$base_img = switch(correct,
            Native = x$ssimg,
            Rigid = x$rig_ssimg
            )
        img = readnii(x$base_img)
        img = window_img(img, 
            window = c(0, 100))

        x$roi_img = switch(correct,
            Native = x$roi,
            Rigid = x$rig_ssroi
            )
        roi = readnii(x$roi_img) > 0.5
        
        # ortho2(img, roi, xyz = xyz)

        pred_img = readnii(img_name) > 0.5

        ###########################
        ### Drop empty image dimensions
        ###########################
        o = dropEmptyImageDimensions(
            img, 
            other.imgs = 
                list(roi = roi, 
                    pred_img = pred_img))
        img = o$outimg
        roi = o$other.imgs$roi
        pred_img = o$other.imgs$pred_img


        # ortho2(img, native.pred, xyz = xyz)
        x = img
        zlim = c(0, 100)
        col = gray(0:64/64)
        xbreaks <- c(min(x, zlim, 
            na.rm = TRUE), 
                  seq(min(zlim, 
                    na.rm = TRUE),
                      max(zlim, 
                        na.rm = TRUE), 
                      length = 
                      length(col) - 1),
                  max(x, zlim, 
                    na.rm = TRUE))


        diff = img
        # false negative
        diff[ roi == 1 & 
            pred_img == 0] = 101
        # false positive
        diff[ roi == 0 & 
            pred_img == 1] = 102
        # true positive
        diff[ roi == 1 & 
            pred_img == 1] = 103

        over_img = niftiarr(img, NA)
        # false negative
        over_img[ roi == 1 & 
            pred_img == 0] = 1
        # false positive
        over_img[ roi == 0 & 
            pred_img == 1] = 2
        # true positive
        over_img[ roi == 1 & 
            pred_img == 1] = 3 
        over_img = cal_img(over_img)

        diff_cols = c("#56B4E9", "#D55E00", 
            "#009E73")
        x = diff
        zlim = c(0, 103)
        col = c(gray(0:64/64),
            diff_cols)
        breaks <- c(min(x, zlim, 
            na.rm = TRUE), 
              seq(min(zlim, na.rm = TRUE),
                  max(zlim, na.rm = TRUE), 
                  length = length(col) - 1),
              max(x, zlim, na.rm = TRUE))
        L = length(breaks)
        breaks[ (L - 2):L] = c(101, 102, 103)
        breaks[ L - 3] = c(100.1)

        plevs = c("False Negative", 
            "False Positive", 
            "True Positive")

        xyz=xyz(roi)
        ### maximal slice
        max_z = which.max(apply(roi, 3, sum))
        zs = which(apply(over_img > 0, 
            3, any))
        med_z = median(zs)


   
        myheight = 3.5
        if (myheight == 7){
            text.y = 20 
        } else {
            text.y = 10
        }
        png(pngname,
            res = 600, units = "in", 
            height= myheight,
            width = 7, type= "cairo")
        double_ortho(img, 
            diff, 
            # don't do alpha blending
            col.y = col,
            xyz = c(xyz[1:2], med_z),
            ybreaks = breaks,
            addlegend = TRUE,
            legend=plevs, 
            leg.col= diff_cols, 
            leg.cex=2, 
            leg.x = -4, leg.y = 45,
            text = sprintf("DSI = %02.2f", 
                run_dice),
            text.cex=1.7, 
            text.y = text.y)
        dev.off()


        png(med_spngname,
            res = 600, 
            units = "in", 
            height= 7,
            width = 7, 
            type= "cairo") 
        overlay(img, over_img,
            plot.type = "single",
            z = med_z,
            col.y = diff_cols)
        dev.off()

        png(max_spngname,
            res = 600, 
            units = "in", 
            height= 7,
            width = 7, 
            type= "cairo") 
        overlay(img, over_img,
            plot.type = "single",
            z = max_z,
            col.y = diff_cols)
        dev.off()        

        # i2 = img[,,c(med_z, max_z)]
        # i2 = copyNIfTIHeader(img = img, 
        #     arr = i2, drop = TRUE)     
        # o2 = over_img[,,c(med_z, max_z)]
        # o2 = copyNIfTIHeader(img = over_img,
        #     arr = o2, drop = TRUE)
        # png(spngname,
        #     res = 600, 
        #     units = "in", 
        #     height=3.5,
        #     width = 7, 
        #     type= "cairo")        
        # overlay(i2, o2,
        #     col.y = diff_cols)   
        # dev.off()
        # spngname = plot_crop(spngname)

        # overlay(img, over_img,
        #     plot.type = "single",
        #     z = med_xyz[3],
        #     col.y = diff_cols)  
                      
        # title(main = "\n\n hey", 
        #     col.main = "white", 
        #     cex.main = 2,
        #     outer = TRUE)
        # legend(x = 0.0, y = 1, 
            # legend = plevs,
        #     fill = diff_cols, 
        #     text.col = "white", bty = "n",
        #     horiz = FALSE,
        #     cex = 1.5)


    }
    print(iimg)
}

