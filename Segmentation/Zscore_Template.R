###################################################################
## This code is for images z-scored compared to a template
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
library(methods)
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

outfile = file.path(outdir, "111_Filenames.Rda")
load(file = outfile)

ttype = "SyN"
interpolator = "Linear"


mean.file = file.path(tempdir, 
  "Mean_Image.nii.gz")
sd.file = file.path(tempdir, 
  "SD_Image.nii.gz")

mean.img = readNIfTI(mean.file)
sd.img = readNIfTI(sd.file)

    
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
if (is.na(iimg)) iimg = 71
vol = fdf$truevol[iimg]


stubfile = function(x, d = NULL, ext = ""){
  b = nii.stub(x, bn=TRUE)
  b = paste0(b, ext)
  file.path(d, b)
}

####################################
## Run both with the Skull Stripped and not skull stripped
####################################

outprefix = tempfile()

x = fdf[iimg,]

ofile = stubfile( x$ss, x$outdir, 
    ext = "_template_zscore.nii.gz" )
# outprefix = stubfile(x$ss, d = x$outdir)

ss.res = ants_regwrite(filename = x$ss, 
    correct = FALSE, 
    outfile = ofile, 
    retimg = TRUE, 
    typeofTransform = ttype,
    template.file = ss.tempfile, 
    interpolator = interpolator,
    remove.warp = FALSE,
    outprefix=outprefix)

z.img = niftiarr(ss.res, (ss.res - mean.img)/ sd.img)
z.img[ is.nan(z.img)] = NA
z.img[ is.infinite(z.img)] = NA

inv.trans = c(
    paste0(outprefix, "0GenericAffine.mat"),
    paste0(outprefix, "1InverseWarp.nii.gz")
    )

z.img[ is.na(z.img)] = 0
ants.ss = antsImageRead(x$ss, 3)
z.img = cal_img(z.img)
z.fname = tempimg(z.img)
ants.z = antsImageRead(z.fname, 3)

zout = antsApplyTransforms(fixed = ants.ss, 
    moving = ants.z, 
    transformlist = inv.trans)

antsImageWrite(zout, ofile)

znative = readNIfTI(ofile, reorient=FALSE)
znative = cal_img(znative)
znative = robust_window(znative)
writeNIfTI(znative, filename = nii.stub(ofile))
ss = readNIfTI(x$ss, reorient=FALSE)

