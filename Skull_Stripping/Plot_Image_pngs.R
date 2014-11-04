###################################################################
## This code is for longitudinal skull stripping volume estimate
##
## Author: John Muschelli
###################################################################
###################################################################
rm(list=ls())
library(plyr)
library(cttools)
library(fslr)
library(matrixStats)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")
rundir = file.path(progdir, "Skull_Stripping")

rerun = FALSE

ids = list.files(basedir, pattern = "^\\d.*\\d$")
iddirs = file.path(basedir, ids)
ssdirs = file.path(iddirs, "Skull_Stripped")
iid <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(iid)) iid = 23
# for (iid in seq_along(ids)){
	id = ids[iid]
	iddir = iddirs[iid]
	outdir = file.path(iddir, "pngs")
	if (!file.exists(outdir)){
		dir.create(outdir)
	}
	imgs = list.files(iddir, pattern = ".nii.gz")
	stubs = nii.stub(imgs)
	imgs = file.path(iddir, imgs)
	pngnames = file.path(outdir, paste0(stubs, ".png"))
	iimg = 1
	for (iimg in seq_along(imgs)){
		pngname = pngnames[iimg]
		if (!file.exists(pngname) | rerun){
			img = readNIfTI(imgs[iimg], reorient=FALSE)
			img = window_img(img, replace = "window")
			dz = dim(img)[3]
			if (dz > 1){
				png(pngname, type="cairo")
				ortho2(img, text = stubs[iimg], 
					text.cex = 1)
				dev.off()
			}
		}
		print(iimg)
	}


	# save(df, file=file.path(outdir, "Skull_Strip_Volumes.Rda"))
# }