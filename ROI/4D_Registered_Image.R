rm(list=ls())
library(oro.nifti)
library(plyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(brainR)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
source(file.path(progdir, "convert_DICOM.R"))
source(file.path(progdir, "fslhd.R"))
basedir = file.path(rootdir, "Registration")


#### baseline NIHSSS ata 
nihss = read.csv(file.path(basedir, "baseline_NIHSS.csv"), 
                 stringsAsFactors=FALSE)
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

outdir = file.path(basedir, "results")
whichdir = "reoriented"

# temp = file.path(tempdir, "scct_unsmooth.nii.gz")
# template = readNIfTI(temp)
# dtemp = dim(template)


cutdown = function(img, outdim= c(181, 217, 181)){
	x = img@.Data 
	x = x[1:outdim[1], 1:outdim[2], 1:outdim[3]]
	img@.Data = x
	# dim(img) = outdim
	img@dim_ = c(3, outdim, 1, 1, 1, 1)
	return(img)
}

#### get all roi images
niis = list.files(pattern='ROI.*.nii', 
	full.names=TRUE, path=basedir, recursive=TRUE)

df = data.frame(dir = dirname(niis), fname = basename(niis),
	stringsAsFactors = FALSE)
df = df[ !grepl("2mm_", df$fname), ]
if (grepl("reoriented", whichdir)) {
	df = df[ grepl("bws", df$fname), ]
}
df = df[ !grepl("affine9", df$fname), ]

df$id = df$dir
df$id = gsub(basedir, "", df$id)
df$id = gsub(".*(\\d\\d\\d-(\\d|)\\d\\d\\d)(.*)", "\\1", df$id)
df$id = as.numeric(gsub("-", "", df$id))

df = df[ df$id == "100318", ]
df$fname = file.path(df$dir, df$fname)


##### get MNI 1mm template
template <- readNIfTI(system.file("MNI152_T1_1mm_brain.nii.gz", 
	package="brainR"), reorient=TRUE)

gm_mask = readNIfTI(file.path(tempdir, "grey_1mm.nii.gz"))
# template[ gm_mask < 0.5] = 0
template=  cutdown(template)
gm_mask =  cutdown(gm_mask)
dtemp <- dim(template)


### 4500 - value that empirically value that presented a brain with gyri
### lower values result in a smoother surface
tt = template
# tt[ tt < 4900 | tt > 6800 ] = 0
orthographic(tt)
# tt[ tt >= 4900 & tt <= 6000 ] = 1
brain <- contour3d(tt, x=1:dtemp[1], y=1:dtemp[2],
z=1:dtemp[3], level = 4500, alpha = 0.2, draw = FALSE, color="black")

files = df$fname
# nii = llply(imgs, readNIfTI, reorient= TRUE, .progress= "text")

### Each visit is a binary mask of hemorrhage in the brain
scene <- list(brain)
## loop through images and thresh
nimgs <- length(files)
cols <- rep("red", nimgs)
iimg = 1
pb = txtProgressBar(min=0, max=1, style=3)
for (iimg in seq(nimgs)) {
	mask <- readNIfTI(files[iimg], reorient=TRUE)
	mask <- drop(mask)
	mask[ is.na(mask) ] = 0
	if (sum(mask > 0) == 0) next;
    setTxtProgressBar(pb, iimg/nimgs)

	### use 0.99 for level of mask - binary
	activation <- contour3d(mask, level = c(0.99), alpha = 1,
	add = TRUE, color=cols[iimg], smooth = 2, draw=FALSE)
	## add these triangles to the list
	scene <- c(scene, list(activation))
}
## make output image names from image names
imgs = basename(files)
fnames <- c("brain.stl", gsub(".nii", ".stl", files, fixed=TRUE))
outfile <-  file.path(outdir, "100-318_4D_outcome.html")

caps = gsub("^bws100-318_", "", imgs)
caps = gsub("(.*)_CT_\\d.*", "\\1", caps)
write4D(scene=scene, fnames=fnames, outfile=outfile,
	captions = c("brain", caps),
	standalone=TRUE, rescale=TRUE, reprint=FALSE, toggle = "radio")
# browseURL(outfile)


