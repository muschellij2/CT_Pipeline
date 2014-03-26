#################################
# Creates map of top p-values from unadusted model
# Author: John Muschelli
#################################
rm(list=ls())
library(oro.nifti)
library(plyr)
library(scales)
library(RColorBrewer)
library(data.table)
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
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

outdir = file.path(basedir, "results")

whichdir = "reoriented"


outfile = file.path(outdir, "Pvalue_Matrix.Rda")
load(file=outfile)

#### load voxel data
outfile = file.path(outdir, "Voxel_Matrix.Rda")
vmat = load(file=outfile )

outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile)

ms = mat[rs > 0, ]
dist.mat = t(ms)
dist.mat = dist.mat %*% ms
A = matrix(diag(dist.mat), ncol=ncol(dist.mat), nrow=nrow(dist.mat))
At = t(A)

dice = 2* dist.mat / (At + A)
stopifnot(all(diag(dice) == 1))

### get top N voxels
y = res[,"X","mod.1"]
rrn = which(rs > ncut)
rrn = rrn[order(y)]



#### reading in template
template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)
dtemp = dim(temp)



open.dev = function(file, type= "cairo", ...){
	get.ext = gsub("(.*)\\.(.*)$", "\\2", file)
	stopifnot( get.ext %in% c("pdf", "bmp", "svg", "png", 
		"jpg", "jpeg", "tiff"))

	## device is jpeg
	if (get.ext == "jpg") get.ext = "jpeg"
	### difff arguments for diff devices
	if (get.ext %in% c("pdf")) {
		do.call(get.ext, list(file=file, ...))
	} else if (get.ext %in% 
		c("bmp", "jpeg", "png", "tiff", "svg")) {
		do.call(get.ext, list(filename=file, type= type,...))
	}
}

view.png = function(fname){
	system(sprintf("display %s", fname))
}

device = "png"

for (nkeep in c(30, 100, 500, 1000)){
	rn = rrn[seq(nkeep)]

	fp = file.path(outdir, 
		paste0("Top_", nkeep, "_pvalues", ".", device))
	open.dev(fp)
	res.p = temp
	res.p[!is.na(res.p)] = NA
	res.p[rn] = 1

	xyz= NULL
	xyz = floor(colMeans(which(res.p > 0, arr.ind = TRUE)))

	# res.p[ rs > ncut ]  = 1-y
	# cols = c("blue", "green", "yellow", "orange", "red", "white")
	mask.overlay(temp, res.p, col.y="red", xyz= xyz)
	dev.off()
}


nkeep = 1000
rn = rrn[seq(nkeep)]

submat = mat[rn,]
wi = colSums(submat)


outfile = file.path(outdir, 
	paste0("Top_", nkeep, "_Pvalues_df.Rda"))
save(submat, rs, rn, wi, nkeep, 
	dist.mat, dice, 
	file=outfile)
