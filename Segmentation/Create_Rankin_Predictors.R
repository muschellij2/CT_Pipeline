###########################################################
## This creates the predictors for the rankin Score
###########################################################
###########################################################
rm(list=ls())
library(plyr)
library(cttools)
library(fslr)
library(ROCR)
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
segdir = file.path(progdir, "Segmentation")


#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )

outfile = file.path(outdir, 
    "111_Filenames_with_volumes_stats.Rda")
load(file = outfile)

fdf$hdr = file.path(fdf$iddir, "Sorted",
    paste0(nii.stub(fdf$img, bn=TRUE), 
        "_Header_Info.Rda"))

iimg <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(iimg)) iimg = 1


fdf$brain_vol = fdf$mean_HU = 
	fdf$median_HU = fdf$sd_HU = NA

for (iimg in seq(nrow(fdf))){

    x = fdf[iimg,] 
    vol = x$truevol

    img = readNIfTI(x$img, reorient = FALSE)
    mask = readNIfTI(x$mask, reorient = FALSE)
    roi = readNIfTI(x$roi, reorient = FALSE)

    vals = img[roi == 1]
    fdf$mean_HU[iimg] = mean(vals)
    fdf$median_HU[iimg] = median(vals)
    fdf$sd_HU[iimg] = sd(vals)


    ####################################
    ## Run both with the Skull Stripped and 
    ## not skull stripped
    ####################################
    x$preddir = x$outdir
    hdrload = load(x$hdr)


    fdf$brain_vol[iimg] = get_roi_vol(mask, 
        dcmtables, mask = mask)$truevol
    print(iimg)
}

ffdf = fdf[, c("id", "mean_HU", "median_HU", 
	"sd_HU", "brain_vol", "truevol")]
outfile = file.path(outdir, 
	"Potential_Rankin_Imaging_Predictors.csv")
write.csv(ffdf, file = outfile, row.names = FALSE)
