#####################################
## Author: John Muschelli
## Date: January 20, 2014
## Purpose: Read in the AAL atlas labels and make R objects
## that can be used later for overlap metrics.
#####################################
#####################################
rm(list=ls())
library(R.matlab)
library(oro.nifti)
library(plyr)
homedir = "/Applications"
basedir = "/Volumes/DATA_LOCAL/Image_Processing/Test_Registration"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  basedir = "/dexter/disk2/smart/stroke_ct/ident/Test_Registration"
}
spmdir = file.path(homedir, "spm8")
aaldir = file.path(spmdir, "toolbox", "aal_for_SPM8")

idir = 1
for (idir in 1:2){
	whichdir = c("reoriented", "FLIRT")[idir]
	# whichdir = "reoriented"
	redir = file.path(basedir, whichdir)
	imgdir <- file.path(redir, "results")

	outfile = file.path(redir, "ROI_Overlap_Measures.Rda")
	file.remove(outfile)
	load(file=file.path(aaldir, "ROIList.Rda"))

	imgs = list.files(path=redir, 
		full.names=TRUE, 
		recursive=FALSE, 
		pattern="^2mm_.*ROI.*")

	imgs = imgs[!grepl("affine9", imgs)]

	iimg = 1;
	x = imgs[iimg]
	# for (iimg in seq_along(imgs)){
	res.list = lapply(imgs, function(x){
		print(x)
		img.fname = x
		img = readNIfTI(img.fname)
		res = get.pct(img.fname, keepall=TRUE)
		cres = collapse.res(res, add.binval=FALSE)
		cres$fname = x
		iimg <<- iimg + 1
		return(cres)
	})

	allres = do.call("rbind", res.list)
	allres$fname = basename(allres$fname)
	allres$fname = gsub("\\.gz$", "", allres$fname)
	allres$fname = gsub("\\.nii$", "", allres$fname)
	allres$fname = gsub("^2mm_", "", allres$fname)
	allres$fname = gsub("affine(9|12)_", "", allres$fname)

	rois$area = rois$name
	allres = merge(allres, rois[, c("area", "col_name")], by="area", 
		all.x=TRUE, sort=FALSE)
	save(allres, res.list, 
		file=outfile)
}
	# cres = cres[order(cres$weighted),]
# }