rm(list=ls())
library(R.matlab)
library(oro.nifti)
library(plyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(knitrBootstrap)
library(xtable)
library(scales)
library(fslr)
library(cttools)
### need cairo for cluster
options(bitmapType='cairo')
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
# source(file.path(progdir, "convert_DICOM.R"))
# source(file.path(progdir, "fslhd.R"))
basedir = file.path(rootdir, "Registration")


#### baseline NIHSSS ata 
nihss = read.csv(file.path(basedir, "baseline_NIHSS.csv"), 
                 stringsAsFactors=FALSE)
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")


outdir = file.path(basedir, "QC")
whichdir = "reoriented"


all.ids = list.dirs(basedir, recursive=FALSE, full.names=FALSE)
all.ids = all.ids[grepl("\\d\\d\\d-(\\d|)\\d\\d\\d", all.ids)]
redir = file.path(basedir, all.ids, whichdir)
redir = redir[file.exists(redir)]
nfiles = sapply(redir, function(x) 
	length(dir(path=x, pattern="ROI*.nii.*")))
redir = redir[nfiles > 0]
## those that have ROIs
all.ids = gsub(paste0(basedir, "/(.*)/", whichdir), "\\1", redir)
# all.ids = all.ids[file.exists(redir)]
uid = as.numeric(gsub("-", "", all.ids))


iid = 1;

for (iid in seq_along(all.ids)){
	

	ID = all.ids[iid]
	iddir = file.path(basedir, ID)
	redir = file.path(iddir, whichdir)

	imgs = list.files(pattern='ROI.*.nii', full.names=TRUE, path=redir)
	imgs = imgs[ !grepl("2mm_", imgs)]
	if (grepl("reoriented", whichdir)) imgs = imgs[ grepl("bws", imgs)]
	imgs = imgs[ !grepl("affine9", imgs)]

	rawimgs = gsub("bws", "w", imgs)
	rawimgs = gsub("^ROI_", "", rawimgs)
	rawimgs = gsub("ROI\\.", ".", rawimgs)


	iimg = 1;
	# for (iimg in seq(nimgs)){
		# if (verbose) print(iimg)

	id = basename(rawimgs[iimg])
	id = nii.stub(id)
	id = gsub("^w", "", id)

	nonreg.imgs = paste0(id, ".nii.gz")
	nonreg.imgs = file.path(iddir, nonreg.imgs)

	rawimgs = rawimgs[iimg]
	file.copy(nonreg.imgs, outdir)
	
	print(iid)
}