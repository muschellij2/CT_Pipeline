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

setwd(aaldir)
x = readMat('ROI_MNI_V4_List.mat')
rois = t(array(x$ROI, dim=c(3, 116)))
rois <- data.frame(rois, stringsAsFactors=FALSE)
names(rois) <- c("label", "short.name", "name")
for (icol in 1:ncol(rois)) {
	rois[, icol] = unlist(rois[, icol])
}
### adding an uncategorized value
rois = rbind(rois, c(0, "UNK", "Uncategorized"))
rois$label = as.numeric(rois$label)

atlas = readNIfTI('ROI_MNI_V4.nii')
at.dim = dim(atlas)
at.pdim = pixdim(atlas)[2:4]
iroi = 1;
ind.list = lapply(rois$label, function(x){
	print(x)
	which(atlas == x)
})
names(ind.list) = rois$name
ind.lengths = sapply(ind.list, length)
all.inds = sort(unique(unlist(ind.list)))


tab.area = function(binimg, keepall) {
	## get overlap of indices
	raw.ind = which(binimg)
	raw.mat = sapply(ind.list, function(x) raw.ind %in% x)
	## cs is sum of indices of overlap
	cs.raw = colSums(raw.mat)
	if (!keepall) cs.raw = cs.raw[cs.raw != 0]
	cs.raw = data.frame(nvox=cs.raw)	
	return(cs.raw)
}
bin.val = 0.95

get.pct = function(img, bin.val = 0.95, rescale=TRUE, keepall=FALSE){
	## you can pass either nifti iamge or filename 
	if (!inherits(img, 'nifti')){
		img = readNIfTI(img)
	}

### make sure same image

	stopifnot(all( at.pdim == pixdim(img)[2:4]))
	stopifnot(all( at.dim == dim(img)))

	## remove NAs and NaN - put to 0
	img[is.na(img)] = 0
	stopifnot(all(img >= 0))
	stopifnot(all(is.finite(img)))
	# img[is.nan(img)] = 0
	###
	img.sum = sum(img)

	## get binary image (any greater than 0)
	rawbin.img = img > 0
	n.vox = sum(rawbin.img)

	### you may want to only call it a true voxel if greater than 
	## bin.val (because of resampling)
	bin.img = img > bin.val
	nbin.vox = sum(bin.img)

	## get overlap
	raw.tab = tab.area(rawbin.img, keepall=keepall)
	bin.tab = tab.area(bin.img, keepall= keepall)

## make sure everything went correctly
	stopifnot(sum(raw.tab$nvox) == n.vox)
	stopifnot(sum(bin.tab$nvox) == nbin.vox)

	img.ind = which(rawbin.img)
	### these voxels have no categorization
	img.mat = sapply(ind.list, function(x) img.ind %in% x)
	img.mat = img.mat[,rownames(raw.tab), drop=FALSE]
	### weighted sum of area
	wsum = apply(img.mat, 2, function(inds){
		grab.ind = img.ind[inds]
		sum(img[grab.ind])
	})

	tol = 1e-5
	diff = abs(sum(wsum) - img.sum)
	stopifnot(diff < tol)
	wsum = data.frame(nvox=wsum)
	if (rescale){
		raw.tab = raw.tab 	/ n.vox
		bin.tab = bin.tab 	/ nbin.vox
		wsum 	= wsum 		/ img.sum		
	}

	return(list(
		bin.val			=	bin.val, 
		weighted.sum	=	wsum, 
		raw.tab			=	raw.tab,
		bin.tab			=	bin.tab,
		n.vox 			=	n.vox,
		nbin.vox		= 	nbin.vox,
		img.sum 		= 	img.sum))
}


collapse.res = function(res, add.binval=FALSE){
	cc = res[c("raw.tab", "bin.tab", "weighted.sum")]
	cc = lapply(cc, function(x){
		x$area = rownames(x)
		rownames(x) = NULL
		return(x)
	})
	c.res = merge(cc$raw.tab, cc$bin.tab, by="area", 
		all=TRUE, suffixes=c(".raw", ".bin"))
	c.res = merge(c.res, cc$weighted.sum, by="area", all=TRUE)
	colnames(c.res) = c("area", "raw", "bin", "weighted")
	nas = sapply(c.res, is.na)
	c.res[nas] = 0
	if (add.binval) c.res$bin.val = res$bin.val
	return(c.res)	
}

roinames = rois$name
roinames = gsub("(.*)_(L|R)$", "\\1", roinames)
roinames = gsub("(.*)_\\d$", "\\1", roinames)
roinames = gsub("(.*)_\\d(|\\d|b)$", "\\1", roinames)
roinames = gsub("(.*)_Crus(1|2)$", "\\1", roinames)
roinames = gsub("Cerebelum", "Cerebellum", roinames)
roinames = gsub("(.*)_(Sup|Inf)$", "\\1", roinames)
roinames = gsub("(.*)_(Mid|Post|Ant)$", "\\1", roinames)
roinames = gsub("(.*)_(Medial|Med|Orb|Tri|Oper)$", "\\1", roinames)
roinames = gsub("(.*)_(Sup|Inf)$", "\\1", roinames)
roinames = gsub("(.*)_(Medial|Med|Orb|Tri|Oper)$", "\\1", roinames)

rois$col_name = roinames

save(rois, ind.list, at.dim, at.pdim,
	ind.lengths, all.inds, 
	get.pct, tab.area, 
	collapse.res,
	file=file.path(aaldir, "ROIList.Rda"))


# img = readNIfTI("/dexter/disk2/smart/stroke_ct/ident/Test_Registration/FLIRT/2mm_affine12_100318_ROI_20010723_0956.nii.gz")



