 
rm(list=ls())
library(R.matlab)
library(oro.nifti)
library(plyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(knitrBootstrap)
library(xtable)
library(brainR)
homedir = "/Applications"
rootdir = "~/CT_Registration"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

outdir = file.path(basedir, "results")
whichdir = "reoriented"

rerun = FALSE

# tempfile = file.path(tempdir, "scct_unsmooth.nii.gz")

tempfile = file.path("~/CT_Registration/clinicaltoolbox/high_res/Skull_Stripped/scct_unsmooth_SS_First_Pass_0.1.nii.gz")

template = readNIfTI(tempfile)

masked = template
# masked[ masked < 0 | masked > 100 ] = 0
dtemp = dim(template)

 brain <- contour3d(masked, x=1:dtemp[1], y=1:dtemp[2],
     z=1:dtemp[3], level = 0.1, alpha = 0.1, draw = FALSE)

# histfile = file.path(outdir, paste0(whichdir, "_Weighted_Sum_Image.nii.gz"))
histfile = file.path("~/CT_Registration", paste0(whichdir, "_Weighted_Sum_Image.nii.gz"))

histimg = readNIfTI(histfile)

x = makeScene(histimg, cutoffs = c(1, 5, 10), alpha= c(0.1, 0.1, 0.1), 
	cols = c("blue", "purple", "red"))

scene = c(list(brain), x)

write4D(scene, outfile = "~/CT_Registration/100_Report_Figure.html", fnames = c("brain.obj", "thresh1.obj", "thresh5.obj", "thresh10.obj"))


