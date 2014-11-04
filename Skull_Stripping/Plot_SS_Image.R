###################################################################
## This code is for longitudinal skull stripping volume estimate
##
## Author: John Muschelli
###################################################################
###################################################################
rm(list=ls())
library(fslr)
library(cttools)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
rerun = FALSE

progdir = file.path(rootdir, "programs")
basedir = file.path(rootdir, "Registration")
datadir = file.path(basedir, "data")
resdir = file.path(basedir, "results")
outdir = file.path(resdir, "SS_Overlays")
if(!file.exists(outdir)){
	dir.create(outdir)
}
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")
rundir = file.path(progdir, "Skull_Stripping")

outfile = file.path(resdir, "Longitudinal_Skull_Strip_Data.Rda")
load(file= outfile)

sdf = split(all.df, all.df$img)

idf <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(idf)) idf = 1133

# for (idf in seq_along(sdf)){
	df = sdf[[idf]]
	img = readNIfTI(df$img[1], reorient=FALSE)
	iimg = 1
	for (iimg in seq(nrow(df))){
		fname = df$fname[iimg]
		pngname = file.path(outdir, 
			paste0(nii.stub(fname, bn = TRUE), '.png'))
		if (!file.exists(pngname) | rerun){

			mask = readNIfTI(fname, reorient=FALSE)

			png(pngname, type="cairo")
				mask.overlay(img, mask, 
					window=c(0, 100))
			dev.off()
		} 
		print(iimg)
	}
	print(idf)
# }