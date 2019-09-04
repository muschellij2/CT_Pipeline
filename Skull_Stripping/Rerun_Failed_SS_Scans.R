###################################################################
## This code is for rerunning failed ss scans, with filling and neck
## removal
##
## Author: John Muschelli
###################################################################
###################################################################
rm(list=ls())
library(plyr)
library(cttools)
library(fslr)
library(extrantsr)
library(methods)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
basedir = file.path(rootdir, "Registration")
datadir = file.path(basedir, "data")
resdir = file.path(basedir, "results")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")
rundir = file.path(progdir, "Skull_Stripping")

########################################
# Read in Failed data
########################################
fail = read.csv(file = file.path(datadir, "Failed_SS_Scans.csv"),
	stringsAsFactors=FALSE)

fail = fail[ fail$Good < 0.25, ]
fail$ssfile = gsub("_Mask", "", fail$ssimg)

fail$hdr = file.path(dirname(fail$img), "Sorted", 
	paste0(nii.stub(fail$img, bn=TRUE), "_Header_Info.Rda"))

cfilters = sapply(fail$hdr, function(hdr){
	load(hdr)
	filt.col = grep("onvolu", colnames(dcmtables), value=TRUE)
	print(filt.col)
	if (length(filt.col) > 0){
		cfilter = unique(dcmtables[, "0018-1210-ConvolutionKernel"])
	} else {
		cfilter = NA
	}
	cfilter
})
get.height = function(file, ...){
	vheight = as.numeric(fslval(file, "pixdim3", ...))
	height = as.numeric(fslval(file, "dim3", ...))
	return(vheight * height)
}

get.vol = function(file, ...){
	vol = fslstats(file, opts = "-l 0 -u 100 -V ", 
		...)
	vol = as.numeric(strsplit(vol, " ")[[1]])
	return(vol)
}

get.sd = function(file, ...){
	vol = fslstats(file, opts = "-s -l 0 -u 100 -S ", 
		...)
	vol = as.numeric(strsplit(vol, " ")[[1]])
	return(vol)
}
# heights = sapply(fail$img, get.height)
# vols = sapply(fail$img, get.vol)
# sds = sapply(fail$img, get.sd)


iimg <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(iimg)) iimg = 22

pngnames = NULL

for (iimg in seq(nrow(fail))){
	print(iimg)
	file = fail[iimg, "img"]

	img = readNIfTI(file, reorient=FALSE)

	ssfile = fail[iimg, "ssfile"]
	ssmaskfile = fail[iimg, "ssimg"]
	ssimg = readNIfTI(ssmaskfile, reorient=FALSE)
	pngname = file.path(dirname(file), 
		"Skull_Stripped", "overlays", 
		paste0(nii.stub(ssmaskfile, bn=TRUE), '.png'))
	png(pngname, type="cairo")
	mask.overlay(img, ssimg, window=c(0, 100))
	dev.off()
	pngnames = c(pngnames, pngname)
}
	# nvoxels =3 



	lthresh = 0
	uthresh = 100

	thresh = niftiarr(img, img * (img > lthresh & img < uthresh))
	## need to do rep.value = 0 because template is like that.
	# noneck = remove_neck(noneck, rep.value=0)
	# reneck = img
	# while(!all((noneck-reneck) == 0)){
		# reneck = noneck
		# noneck = remove_neck(reneck, rep.value=0)
	neck_mask = remove_neck(thresh, rep.value=0, 
		template.file = 
			system.file("scct_unsmooth_Skull_Stripped.nii.gz", 
			package = "ichseg"),
		ret_mask = TRUE	
		)
	noneck = mask_img(img, neck_mask)
	# apply(noneck, 3, function(x) all(x == 0))
	ortho2(noneck, window=c(0, 100))
	# }

	# ofile = tempfile()
	ofile = ssfile
	ss = CT_Skull_Strip(noneck, outfile = ofile, retimg = TRUE)
	ssmask = readNIfTI(paste0(nii.stub(ofile), "_Mask"), 
		reorient = FALSE)
	mask.overlay(img, ssmask, window=c(0, 100))

	old.cog = cog(img > 0 & img < 100, ceil=TRUE)
	cog = cog(ssmask, ceil=TRUE)
	ss = CT_Skull_Strip(noneck, outfile = ofile, retimg = TRUE, 
		opts = paste("-f 0.01 -v -w 2", 
			paste(c("-c", cog), collapse=" ")))
	ssmask = readNIfTI(paste0(nii.stub(ofile), "_Mask"), 
		reorient = FALSE)
	mask.overlay(img, ssmask, window=c(0, 100))


	ss_fill = dil_ero(ssmask, retimg= TRUE, nvoxels =5)
	mask.overlay(img, ss_fill, window=c(0, 100))

	ss_img = mask_img(img, ss_fill)
	ofile = nii.stub(ofile)
	writeNIfTI(ss_img, file = ofile)

	omaskfile = nii.stub(ssmaskfile)
	writeNIfTI(ss_fill, file = omaskfile)

# }



# # mask.overlay(img, ssmask, window=c(0, 100))

# ss_fill = dil_ero(ssmask, retimg= TRUE, nvoxels =5)
# # ss_fill2 = dil_ero(ssmask, retimg= TRUE, nvoxels =5, 
# #shift = FALSE)

# ss_fill2 = dil_ero_mmand(ssmask, retimg= TRUE, nvoxels =5)
# all.equal(ss_fill, ss_fill2)



