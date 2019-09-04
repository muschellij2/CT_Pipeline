rm(list=ls())
library(fslr)
library(WhiteStripe)
library(plyr)
library(stringr)
library(data.table)
library(cttools)
library(extrantsr)
library(methods)
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

iddir = file.path(basedir, "179-395")
setwd(iddir)

outdir = file.path(iddir, "catheters")
if (!file.exists(outdir)){
	dir.create(outdir)
}

imgs = list.files(path = iddir, full.names=TRUE, pattern='.nii.gz')
ss.tempfile = file.path(tempdir, "Skull_Stripped",
    "scct_unsmooth_SS_0.01.nii.gz")
ss.tempmask = file.path(tempdir, "Skull_Stripped",
    "scct_unsmooth_SS_0.01_Mask.nii.gz")

iimg = 1
uthresh = 90
mask_prob = 0.99
boxsize = 5
rerun = TRUE


for (iimg in seq_along(imgs)){

	img.f = imgs[iimg]

	outfile = file.path(outdir,
		paste0(nii.stub(img.f, bn=TRUE), "_cath.nii.gz"))

	if (!file.exists(outfile) | rerun){

		img = readNIfTI(img.f, reorient=FALSE)

		ssfile = file.path(iddir, 
			'Skull_Stripped',
			paste0(nii.stub(img.f, bn=TRUE), "_SS_0.01.nii.gz"))

		ssmaskfile = file.path(iddir, 
			'Skull_Stripped',
			paste0(nii.stub(img.f, bn=TRUE), 
				"_SS_0.01_Mask.nii.gz"))

		ssmask = readNIfTI(ssmaskfile, reorient=FALSE)

		ofile_ss = paste0(tempfile(), '.nii.gz')
		ofile_ssmask = paste0(tempfile(), '.nii.gz')

		oss = ants_regwrite(ss.tempfile, 
				outfile = ofile_ss, 
		       template.file = ssfile, interpolator = "Linear",
		       other.files = ss.tempmask, 
		       other.outfiles = ofile_ssmask,
		       typeofTransform = "SyN",
		       remove.warp = TRUE)
		mask = readNIfTI(ofile_ssmask, reorient=FALSE)

		mask = cal_img(mask > mask_prob)

		ero_mask = fslerode(file=mask, 
			kopts = paste0("-kernel boxv ", boxsize), 
			reorient=FALSE, retimg = TRUE)

	## Of the SS mask, those that are in the template mask but not
	# the
		## SS mask

	# mask.overlay(img, ero_mask, col.y="red", window=c(0, 100))

		diff_mask = niftiarr(ssmask, (ero_mask - ssmask) > 0)

		# masked = mask_img(img, diff_mask)
		masked = mask_img(img, ero_mask)

		cath = cal_img(masked > uthresh)
		ofile = paste0(tempfile(), '.nii')
		spm_bwlabel(cath, outfile = ofile, k = 1500)
		spm_bwlabel(ofile, outfile = ofile, topN = 1)
		ocath = readNIfTI(ofile, reorient=FALSE)
	} else{
		fslany = function(...){
			as.numeric(fslstats(..., opts = "-M")) >0
		}
		print(fslany(outfile))
	}

	# mask.overlay(img, ocath, col.y="red", window=c(0, 100))
}

# maskfile = file.path(iddir, 'Skull_Stripped", 
# 	"179-395_20110113_1541_CT_80364_CT_AXIAL_SS_0.01_Mask.nii.gz')

# mask = readNIfTI(maskfile, reorient=FALSE)
# outfile = paste0(tempfile(), '.nii')
# imfill(infile = maskfile, outfile = outfile, margin=3)
# out = readNIfTI(outfile, reorient=FALSE)
# ortho2(out)

# slice = out[,,15]

# masked = mask_img(img, out)
# cath = cal_img(masked > 100)

# ofile = paste0(tempfile(), '.nii')
# spm_bwlabel(cath, outfile = ofile, topN = 1)
# ocath = readNIfTI(ofile, reorient=FALSE)

# mask.overlay(img, ocath, col.y="red", window=c(0, 100))


# mask_fill = function(mask){

# 	ind = which(mask > 0, arr.ind=TRUE)
# 	cog = colMeans(ind)
# 	dist = sqrt(rowSums(t(t(ind) - cog)^2))
# 	d = t(t(ind) - cog)
# 	theta = acos(d[, "dim3"] / dist)
# 	phi = atan(d[, "dim2"]/d[, "dim1"])
# 	arr = array(0, dim = dim(mask))
# 	arr[ind] = dist


# 	t.quants = quantile(theta, probs = seq(0, 1, by=0.01))
# 	cut.theta = cut(theta, breaks=t.quants, include.lowest= TRUE)

# 	p.quants = quantile(phi, probs = seq(0, 1, by=0.01))
# 	cut.phi = cut(phi, breaks=p.quants, include.lowest= TRUE)	

# 	df = data.frame(r = dist, 
# 		theta = theta, t.quant = cut.theta, 
# 		phi = phi, p.quant = cut.phi)

# 	maxes = ddply(df, .(t.quant, p.quant), summarise, 
# 		maxr = max(r))
# 	nums = strsplit(as.character(maxes$t.quant), ",") 
# 	nums= lapply(nums, function(x){
# 		x = gsub("[", "", x, fixed=TRUE)
# 		x = gsub("]", "", x, fixed=TRUE)
# 		x = gsub("(", "", x, fixed=TRUE)
# 		x = gsub(")", "", x, fixed=TRUE)
# 		x = str_trim(x)
# 		x = as.numeric(x)
# 		(x[1] + x[2]) /2
# 	})
# 	maxes$t.num = nums

# 	nums = strsplit(as.character(maxes$p.quant), ",") 
# 	nums= lapply(nums, function(x){
# 		x = gsub("[", "", x, fixed=TRUE)
# 		x = gsub("]", "", x, fixed=TRUE)
# 		x = gsub("(", "", x, fixed=TRUE)
# 		x = gsub(")", "", x, fixed=TRUE)
# 		x = str_trim(x)
# 		x = as.numeric(x)
# 		(x[1] + x[2]) /2
# 	})
# 	maxes$p.num = nums	

# 	maxes$fac = paste0(maxes$t.quant, "_", 
# 		maxes$p.quant)
# 	maxes$t.quant = maxes$p.quant = NULL	

# 	maxes = data.table(maxes)
# 	setkey(maxes, "fac")
# 	# plot(theta, dist, pch=".")
# 	# lines(maxes$num, maxes$maxr, type="l")

# 	all.ind = expand.grid(dim1 = 1:nrow(mask), 
# 		dim2= 1:ncol(mask), 
# 		dim3 = 1:nsli(mask))
# 	all.r = sqrt(rowSums(t(t(all.ind) - cog)^2))
# 	all.d = t(t(all.ind) - cog)

# 	all.theta = acos(all.d[, "dim3"] / all.r)
# 	all.phi = atan(all.d[, "dim2"]/all.d[, "dim1"])
# 	all.tcut = cut(all.theta, breaks=t.quants, 
# 		include.lowest= TRUE)
# 	all.pcut = cut(all.phi, breaks=p.quants, 
# 		include.lowest= TRUE)	

# 	all.df = data.frame(t.quant = all.tcut, 
# 		p.quant = all.pcut, 
# 		r = all.r, all.ind)
# 	all.df$fac = paste0(all.df$t.quant, "_", 
# 		all.df$p.quant)
# 	all.df$t.quant = all.df$p.quant = NULL

# 	all.df = data.table(all.df)
# 	setkey(all.df, "fac")

# 	all.df = all.df[maxes, ]

# 	all.df[, keep := r <= maxr]
	
# 	keep.ind = as.matrix(all.df[ which(keep > 0), 
# 		list(dim1, dim2, dim3)])

# 	newmask = array(0, dim = dim(mask))
# 	newmask[keep.ind] = 1
# 	newmask =  niftiarr(mask, newmask)
# }

# # ortho2(out, xyz=c(256, 256, 15))



# slice_fill = function(slice, myprob = 1){

# 	ind = which(slice > 0, arr.ind=TRUE)
# 	cog = colMeans(ind)
# 	dist = sqrt(rowSums(t(t(ind) - cog)^2))
# 	d = t(t(ind) - cog)
# 	theta = atan2(d[, 2], d[,1])
# 	arr = array(0, dim = dim(slice))
# 	arr[ind] = dist


# 	quants = quantile(theta, probs = seq(0, 1, by=0.01))

# 	cut.theta = cut(theta, breaks=quants, include.lowest= TRUE)

# 	df = data.frame(r = dist, theta = theta, quant = cut.theta)

# 	maxes = ddply(df, .(quant), summarise, 
# 		maxr = quantile(r, probs = myprob))
# 	nums = strsplit(as.character(maxes$quant), ",") 
# 	nums= lapply(nums, function(x){
# 		x = gsub("[", "", x, fixed=TRUE)
# 		x = gsub("]", "", x, fixed=TRUE)
# 		x = gsub("(", "", x, fixed=TRUE)
# 		x = gsub(")", "", x, fixed=TRUE)
# 		x = str_trim(x)
# 		x = as.numeric(x)
# 		(x[1] + x[2]) /2
# 	})
# 	maxes$num = nums

# 	plot(theta, dist, pch=".")
# 	lines(maxes$num, maxes$maxr, type="l")

# 	all.ind = expand.grid(x = 1:nrow(slice), y= 1:ncol(slice))
# 	all.r = sqrt(rowSums(t(t(all.ind) - cog)^2))
# 	all.d = t(t(all.ind) - cog)
# 	all.theta = atan2(all.d[, 2], all.d[,1])

# 	all.cut = cut(all.theta, breaks=quants, include.lowest= TRUE)

# 	all.df = data.frame(quant = all.cut, r = all.r, all.ind)
# 	all.df = merge(all.df, maxes, by="quant")
# 	all.df$keep = all.df$r <= all.df$maxr
# 	keep.ind = as.matrix(all.df[all.df$keep, c("x", "y")])

# 	newslice = array(0, dim = dim(slice))
# 	newslice[keep.ind] = 1
# 	newslice
# }


# # tofill = which((out - mask) > 0, arr.ind=TRUE)
# # unique(tofill[, "dim3"])

# for (i in seq(nsli(out))){
# 	newslice = slice_fill(mask[,,i])
# 	out@.Data[,,i] = newslice
# 	print(i)
# }

# img.hist = hist(theta, breaks=200)
# y = img.hist$counts
# x = img.hist$mids

# smoothed = smooth_hist(x, y)