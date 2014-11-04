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
resdir = file.path(basedir, "results", "pngs")
if(!file.exists(resdir)){
	dir.create(resdir)
}
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")
rundir = file.path(progdir, "Skull_Stripping")

rerun = FALSE

ids = list.files(basedir, pattern = "^\\d.*\\d$")
iddirs = file.path(basedir, ids)
ssdirs = file.path(iddirs, "Skull_Stripped")
iid <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(iid)) iid = 23
allpngs = NULL
for (iid in seq_along(ids)){
	iddir = iddirs[iid]
	outdir = file.path(iddir, "pngs")
	imgs = list.files(iddir, pattern = ".nii.gz")
	stubs = nii.stub(imgs)

	pngnames = file.path(outdir, paste0(stubs, ".png"))
	pngnames = pngnames[file.exists(pngnames)]
	# save(df, file=file.path(outdir, "Skull_Strip_Volumes.Rda"))
	allpngs = c(allpngs, pngnames)
}
save(allpngs, file=file.path(resdir, "All_PNG_Paths.Rda"))

file.copy(allpngs, to = resdir)