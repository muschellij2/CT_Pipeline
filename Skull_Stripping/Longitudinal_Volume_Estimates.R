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
if (is.na(iid)) iid = 5

# for (iid in seq_along(ids)){
	id = ids[iid]
	iddir = iddirs[iid]
	ssdir = ssdirs[iid]
	outdir = file.path(iddir, "results")
	if (!file.exists(outdir)){
		dir.create(outdir)
	}
	outfile = file.path(outdir, "Skull_Strip_Volumes.Rda")
	if (!file.exists(outfile) | rerun){
		masks = list.files(ssdir, pattern = "_Mask.nii.gz")
		masks = masks[!grepl("nopresmooth|refill", masks)]
		int = gsub(".*SS_(.*)_Mask.*", "\\1", masks)
		stub = gsub("(.*)_SS_(.*)_Mask.*", "\\1", masks)
		hdr = file.path(iddir, "Sorted", 
			paste0(stub, "_Header_Info.Rda"))
		stopifnot(all(file.exists(hdr)))

		df = data.frame(matrix(NA, nrow=length(masks), 
			ncol= 6))
		colnames(df) = c("fname", "hdr", "truevol", "thickvol", 
			"zvol", "varslice")
		df$fname = file.path(ssdir, masks)
		df$hdr = hdr
		iimg =1 
		for (iimg in seq(nrow(df))){
			load(df$hdr[iimg])
			img = readNIfTI(df$fname[iimg], reorient = FALSE)
			res = get_roi_vol(img = img, dcmtables= dcmtables)
			df[iimg, "truevol"] = res$truevol
			df[iimg, "thickvol"] = res$thickvol
			df[iimg, "zvol"] = res$zvol
			df[iimg, "varslice"] = res$varslice
			print(iimg)
		}
		save(df, file=outfile)
	}
# }