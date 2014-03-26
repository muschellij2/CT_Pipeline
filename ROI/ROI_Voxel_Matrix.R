###############################
# make matrix of p-values for the group
###############################
rm(list=ls())
library(oro.nifti)
library(plyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(knitrBootstrap)
library(xtable)
library(scales)
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

# whichdir = "FLIRT"
rerun = FALSE

template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)
dtemp = dim(temp)



#### get all roi images
imgs = list.files(pattern='ROI.*.nii', 
	full.names=TRUE, path=basedir, recursive=TRUE)

df = data.frame(dir = dirname(imgs), fname = basename(imgs),
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

### grab only the first image
df = ddply(df, .(id), function(x){
	x[1,]
})

nimgs = nrow(df)
mat = matrix(FALSE, nrow=prod(dtemp), ncol=nimgs)
read.img = function(fname){
	img = readNIfTI(fname)
	img[is.nan(img) | is.na(img)] = 0
	# res = which(img > 0)
	# return(res)
	return(c(img > 0))
}
iimg = 1
pb = txtProgressBar(min=0, max=1, style=3)
for (iimg in seq(nimgs)){
	dd = df[iimg, , drop=FALSE]
	mat[, iimg] = read.img(file.path(dd$dir, dd$fname))
	setTxtProgressBar(pb, iimg/nimgs)
}

rs = rowSums(mat)

outfile = file.path(outdir, "Voxel_Matrix.Rda")
save(mat, rs, df, file=outfile )