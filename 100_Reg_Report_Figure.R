 
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
library(scales)
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

# tempfile = file.path("~/CT_Registration/clinicaltoolbox/high_res/Skull_Stripped/scct_unsmooth_SS_First_Pass_0.1.nii.gz")
tempfile = file.path("~/CT_Registration/clinicaltoolbox/high_res//betsct1_unsmooth.nii.gz")

template = readNIfTI(tempfile)

masked = template
# masked[ masked < 0 | masked > 100 ] = 0
dtemp = dim(template)

brain <- contour3d(masked, x=1:dtemp[1], y=1:dtemp[2],
     z=1:dtemp[3], level = 200, alpha = 0.37, draw = FALSE,
     color="black")

# histfile = file.path(outdir, paste0(whichdir, "_Weighted_Sum_Image.nii.gz"))
histfile = file.path("~/CT_Registration/CT_Pipeline/", 
                     paste0(whichdir, "_Binary_Sum_Image.nii.gz"))

histimg = readNIfTI(histfile)
maxn = ceiling(max(histimg)*111)

img_cut = function(img, breaks, ...){
  cuts = cut(img, breaks=breaks, ...)
  # cuts = factor(cuts, levels)
  levs = levels(cuts)
  cuts = as.numeric(cuts)
  # res.p[ rs > ncut ] = cuts
  img@.Data = array(cuts, dim=dim(img))
  return(list(img=img, levs=levs))
}

npts = 111
breaks = seq(0, .45, by=.05)
clist = img_cut(histimg, breaks=breaks, include.lowest=FALSE)
cimg = clist$img
levs = clist$levs
col.cut = div_gradient_pal(low="blue", 
                                 mid="red", 
                                 high="yellow")(
                                   seq(0, 1, length=length(levs))
                                 )
alpha = c(.5, rep(1, length=length(levs)-1))
cimg[is.na(cimg)] = 0
x = makeScene(cimg, cutoffs = seq_along(levs) - .01, alpha=alpha , 
	cols = col.cut)

# x = makeScene(histimg, cutoffs = breaks, alpha=alpha , 
#               cols = c("white", col.cut))

scene = c(list(brain), x)

write4D(scene, 
        outfile = "~/CT_Registration/programs/100_Report_Figure.html", 
        fnames = c("brain.stl", 
                   paste0("thresh", seq_along(levs), ".stl")),
        visible = c(TRUE, rep(FALSE, length(levs))),
        captions = c("brain", paste0("Proportion: ", levs)), 
                     reprint = FALSE, rescale=TRUE, 
        xtkgui = FALSE, toggle="slider")


