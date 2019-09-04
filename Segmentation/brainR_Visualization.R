###################################################
## Registration of all scans to first scan
##
## Author: John Muschelli
## Last updated: May 20, 2014
################################################
##################################################
rm(list=ls())
library(fslr)
library(brainR)
library(getopt)
library(misc3d)
library(cttools)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  ident_rootdir = "/dexter/disk2/smart/stroke_ct/ident"
  # rootdir = "/dexter/disk2/smart/stroke_ct/deident"
  rootdir = ident_rootdir
}
progdir = file.path(rootdir, "programs")
# basedir = file.path(rootdir, "Longitudinal")
# roidir = file.path(rootdir, "Longitudinal_ROI")
basedir = file.path(rootdir, "Registration")
roidir = file.path(rootdir, "ROI_data")
tempdir = file.path(ident_rootdir, "Template")
atlasdir = file.path(ident_rootdir, "atlases")


ids = list.dirs(basedir, recursive=FALSE, 
  full.names=FALSE)
ids = basename(ids)
ids = grep("\\d\\d\\d-(\\d|)\\d\\d\\d", ids, value=TRUE)
length(ids)

ss = TRUE
ttype = "Affine"

iid <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(iid)) iid <- 10

id = ids[iid]
iddir = file.path(basedir, id)
id_roidir = file.path(roidir, id)

#################
# Grab all the nifti files for the skull stripped images
#################
niis = list.files(path = file.path(iddir, 
  ifelse(ss, "Skull_Stripped", "")), 
pattern = paste0(
  ifelse(ss, "SS_0.01", ""), 
  '[.]nii[.]gz'), 
recursive=FALSE,
full.names=TRUE )

basenii = niis[1]
niis = niis[-1]

rois = file.path(id_roidir, 
  paste0(
    gsub(
      ifelse(ss, "_SS_0.01$", "$"),
      "", 
      nii.stub(niis, bn = TRUE)
      ),
    "ROI.nii.gz")
  )

df = data.frame(nii = niis, roi = rois, 
  basenii = basenii,
  stringsAsFactors = FALSE)
df = df[ file.exists(df$roi), , drop = FALSE]

###########################
# Create output directory
###########################
outdir = file.path(iddir, "Coregistered")
if (!file.exists(outdir)){
  dir.create(outdir)
}
df$outdir = outdir
addon = ""
addon = paste0("_", ttype)
df$outfile = file.path(df$outdir, 
  paste0(nii.stub(df$nii, bn=TRUE), 
    "_Coregistered", addon))
df$roi.outfile = file.path(df$outdir, 
  paste0(nii.stub(df$roi, bn=TRUE), 
    "_Coregistered", addon))

ifile = 1

img = readnii(df$basenii[1] )
vdim = voxdim(img)
dimg = mapply(function(x, y){
  round(seq(from=1, to=x)* y)
  }, dim(img), vdim)

mask = img > 0
dmask = fsldilate(mask)
diff_mask = dmask - mask

cc = spm_bwlabel(diff_mask, topN = 1, retimg = TRUE)


brain <- contour3d(cc, x=dimg[[1]], y=dimg[[2]],
  z=dimg[[3]], 
  level = 0.99, alpha = 0.6, smooth = TRUE,
  color = "gray", draw = FALSE)

files = df$roi.outfile
imgs = lapply(files, readnii)
imgs = lapply(imgs, function(x) x > 0.5)
sums = sapply(imgs, sum)
imgs = imgs[sums > 0]
files = files[ sums > 0]

scene <- list(brain)
nimgs <- length(imgs)
cols <- rainbow(nimgs)
for (iimg in 1:nimgs) {
  roi <- imgs[[iimg]]
  #   ss <- strsplit(iimg, "\\.")
  #   stub <- sapply(ss, function(x) x[1])
  
  # mask <- readNIfTI(img_num, reorient=FALSE)
  # if (length(dim(mask)) > 3) mask <- mask[,,,1]
  #mask <- mask[dmask[1]:1, dmask[2]:1,]
  #   mask <- mask[, dmask[2]:1,]
  ## need to reverse the brain to match template
  
  
  ### this would be the ``activation'' or 
  ###surface you want to render 
  ### here just taking the upper WM from 
  ##the template image
  activation <- contour3d(roi, x=dimg[[1]], y=dimg[[2]],
    z=dimg[[3]], level = c(0.99), 
    alpha = 1, 
    add = TRUE, color=cols[iimg], 
    draw=FALSE)
  
  scene <- c(scene, list(activation))
  print(iimg)  
  #   drawScene.rgl(activation)
  
}

fnames <- c("brain.stl", 
  paste0(nii.stub(files, bn = TRUE), ".stl"))
outfile <- file.path(outdir,  "index_4D_stl.html")
write4D(scene=scene, fnames=fnames, 
  outfile=outfile, standalone=TRUE)


