rm(list=ls())
library(fslr)
library(WhiteStripe)
library(plyr)
library(stringr)
library(data.table)
library(cttools)
library(extrantsr)
library(methods)
library(mmand)
library(matrixStats)
### use big protocol for cath
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
basedir = file.path(rootdir, "Registration")
roidir = file.path(rootdir, "ROI_data")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

ids = list.dirs(basedir, recursive=FALSE, full.names=FALSE)
ids = basename(ids)
ids = grep("\\d\\d\\d-(\\d|)\\d\\d\\d", ids, value=TRUE)
length(ids)

iid <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if (is.na(iid)) iid <- 2
# all.imgs = NULL
# for (iid in seq_along(ids)){

	id = ids[iid]
	print(id)
	iddir = file.path(basedir, id)
	# setwd(iddir)

	outdir = file.path(iddir, "catheters")
	if (!file.exists(outdir)){
		dir.create(outdir)
	}

	imgs = list.files(path = iddir, full.names=TRUE, 
		pattern='.nii.gz')

	### need to delete CTAs
	imgs = imgs[ !grepl("CTA", imgs)]

# 	all.imgs = c(all.imgs, imgs)
# 	print(imgs)
# }
ss.tempfile = file.path(tempdir, "Skull_Stripped",
    "scct_unsmooth_SS_0.01.nii.gz")
ss.tempmask = file.path(tempdir, "Skull_Stripped",
    "scct_unsmooth_SS_0.01_Mask.nii.gz")

# iid = 80 and iimg = 4
# iid = 16, iimg = 3
# iid = 16, iimg = 4 has 2 cath
iimg = 10
uthresh = 100
mask_prob = 0.99
boxsize = 7
rerun = TRUE
size = 1000

for (iimg in seq_along(imgs)){

	img.f = imgs[iimg]
	print(img.f)
	outfile = file.path(outdir,
		paste0(nii.stub(img.f, bn=TRUE), 
			"_cath.nii.gz"))

	if (!file.exists(outfile) | rerun){

		img = readNIfTI(img.f, reorient=FALSE)

		ssmaskfile = file.path(iddir, 
			'Skull_Stripped',
			paste0(nii.stub(img.f, bn=TRUE), 
				"_SS_0.01_Mask.nii.gz"))

		ssmask = readNIfTI(ssmaskfile, reorient=FALSE)
		resolution = prod(pixdim(ssmask)[2:4]) / 1000

		z = pixdim(ssmask)[4]

		##### size is 1cc - about 1000 voxels
		size = floor(1/resolution)

		# dil = fslfill2(file = ssmaskfile, 
		# 	kopts = paste0("-kernel box ", boxsize),
		# 	retimg = TRUE)
		nvoxels = 3
		# dil_ero = dil_ero_mmand(ssmask, 
		# 	nvoxels = nvoxels, 
		# 	retimg=TRUE)

		dil_ero = dil_ero(ssmask, nvoxels =nvoxels, 
			retimg=TRUE)
		# box of width 5mm means (voxels * 2 + 1) * dim = 5mm
		# or (5-1)/2 = 2 mm
		# nvoxels = pmax(
		# 	round((boxsize/pixdim(ssmask)[2:4] - 1)/2), 3)
		# dil_ero2 = mean_image(1 - dil_ero, 
		# 	nvoxels = nvoxels) > .Machine$double.eps^0.5
		# dil_ero2 = niftiarr(ssmask, 1 - dil_ero2)

		len = 3
		pdim = pixdim(ssmask)[2:4]
		nz = ceiling(len/pdim)
		nvoxels = nz
		# nvoxels = c(3, 3, nz)

		dilate = mean_image(dil_ero, 
			nvoxels = 1)
		dilate = niftiarr(img, dilate >= 
			(1- .Machine$double.eps^0.5)
			)
		surface = niftiarr(img, dil_ero - dilate)
		surf.ind = which(surface > 0, arr.ind=TRUE)
		surf.mm = t(t(surf.ind) * pdim)
		sq.surf.mm = rowSums(surf.mm^2)
		
		# dil_ero2 = mean_image(dil_ero, 
		# 	nvoxels = nvoxels)
		# ## can't have 10 or more neighbors being 0
		# # 10/(7*7*3)
		# dil_ero2 = niftiarr(ssmask, dil_ero2 >= 
		# 	(1- .Machine$double.eps^0.5)
		# 	)

		len = 1.5
		nz = ceiling(len/pixdim(ssmask)[4])
		nvoxels = c(3, 3, nz)
		dil_ero2 = mean_image(1 - dil_ero, 
			nvoxels = nvoxels)
		dil_ero2 = 1 - dil_ero2
		dil_ero2 = niftiarr(dil_ero, dil_ero2) >= 
			(1- .Machine$double.eps^0.5)

		# box = shapeKernel(width = nvoxels, type="box")
		# dil_ero2 = erode(ssmask, kernel = box)
		# dil_ero2 = niftiarr(ssmask, dil_ero2)

		############### Using FSL
		# kopts = ""
		# if (z < 3){
		# 	## 15 mm does well in z-plane (for 5mm thickness)
		# 	zdim = floor(15/z)
		# 	zdim = 2 *round((zdim+1)/2)-1;
		# 	kopts = paste0("-kernel boxv 3x3x", zdim)
		# }

		# dil = fslfill2(file = ssmaskfile, 
		# 	kopts = kopts,
		# 	retimg = TRUE)

		# ## when from top of head (C trajectory - 
		# need to do within)
		# kopts = "-kernel boxv 7x7x1"
		# dil = fslfill2(file = dil, 
		# 	kopts = kopts,
		# 	retimg = TRUE)		
		# ##### Dilation

		# dil = fslmaths(file = ssmaskfile, 
		# 	opts = "-mul -1 -add 1 -ero -mul -1 -add 1", 
		# 	retimg=TRUE)

		# ###############################
		# # Removing top and bottom slices  should be no cath
		# ################################
		# dil_ero2 = fslerode(file=dil, 
		# 	kopts = paste0("-kernel box ", 
		# 		boxsize, "x", boxsize, "x1"),
		# 	retimg = TRUE)
		# dil_ero2 = fslerode(file=dil2, 
		# 	kopts = paste0("-kernel boxv ", boxsize), 
		# 	reorient=FALSE, retimg = TRUE)


		masked = mask_img(img, dil_ero2)

		cath = cal_img(masked > uthresh)

		########################################
		# Calculating the distance to surface
		########################################
		cath.ind = which(cath > 0, arr.ind=TRUE)
		cath.mm = t(t(cath.ind) * pdim)
		sq.cath.mm = rowSums(cath.mm^2)

		buffsize = 500
		N = nrow(cath.ind)
		nfolds = ceiling(N/buffsize)
		folds = rep(1:nfolds, length.out=N)
		dists = rep(NA, length = N)
		ifold = 1

		for (ifold in seq(nfolds)){
			print(ifold)
			ind = which(folds == ifold)
			cc = cath.mm[ind,]
			sq = sq.cath.mm[ind]
			xy = -2 * tcrossprod(cc, surf.mm)
			d = xy + sq 
			d = t(t(d) + sq.surf.mm)
			d = rowMins(d)
			d[ abs(d) <= sqrt(.Machine$double.eps)] = 0
			dists[ind] = sqrt(d)
		}
		# mask.overlay(masked, cath, col.y="red", window=c(0, 100))


		cath@.Data[ cath.ind[ dists ==0, ]] = 0

		cath.ind = cbind(cath.ind, dist = dists)
		
		ofile = paste0(tempfile(), '.nii')
		spm_bwlabel(cath, outfile = ofile, k = size, binary=FALSE)
		# spm_bwlabel(ofile, outfile = ofile, topN = 1)
		ocath = readNIfTI(ofile, reorient=FALSE)
		ocath = cal_img(ocath)
		
		cath.ind2 = cath.ind[ cath.ind[, "dist"] > 0, ]
		cath.ind2= cbind(cath.ind2, 
			clust=ocath[cath.ind2[, c("dim1", "dim2", "dim3")]])
		bad.ind = cath.ind2[cath.ind2[, "clust"] == 0, 
			c("dim1", "dim2", "dim3")]
		cath@.Data[ bad.ind] = 0

		cath.ind2 = cath.ind2[cath.ind2[, "clust"] > 0,]
		cath.ind2 = data.frame(cath.ind2)
		mean.dist = ddply(cath.ind2, .(clust), summarise,
			mean_dist = mean(dist))

		xyz= cog(ocath, ceil=TRUE)
		mask.overlay(masked, ocath, col.y="red", window=c(0, 100),
			xyz=xyz)



		ofile = paste0(tempfile(), '.nii')
		spm_bwlabel(cath, outfile = ofile, k = size, binary=TRUE)
		# spm_bwlabel(ofile, outfile = ofile, topN = 1)
		ocath = readNIfTI(ofile, reorient=FALSE)
		ocath = cal_img(ocath)
		writeNIfTI(ocath, file=nii.stub(outfile))
		xyz= cog(cath, ceil=TRUE)
		mask.overlay(masked, ocath, col.y="red", window=c(0, 100),
			xyz=xyz)

	} else {
		fslany = function(...){
			as.numeric(fslstats(..., opts = "-M", verbose=FALSE)) >0
		}
		print(fslany(outfile))
	}

	# mask.overlay(img, ocath, col.y="red", window=c(0, 100))
}
