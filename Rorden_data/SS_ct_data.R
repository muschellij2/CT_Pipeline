rm(list=ls())
library(cttools)
library(fslr)
library(R.utils)
homedir = "~/"
rootdir = "~/CT_Registration/"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
basedir = file.path(rootdir, "Registration")

homedir = file.path(rootdir, "Rorden_data")
regdir = file.path(homedir, "CT_lesions")


# scandir = file.path(homedir, "CT_lesions")
scandir = file.path(homedir, "ct_scans")
# setwd(homedir)
outdir = tempdir()


files = list.files(scandir, pattern = "^rx.*.nii", full.names=TRUE)
files = files[ !grepl("_SS_", files)]
ifile = 30


for (ifile in seq_along(files)){
	print(ifile)
	myfile = files[ifile]
	img = clinical_c2h(myfile)
	ssfile = paste0(nii.stub(myfile), "_SS_0.01")
	ss = CT_Skull_Strip(img, 
		retimg=FALSE, 
		lthresh = -50,
		uthresh = 100,
		sigma = 3,
		outfile = ssfile)
	# ssfile = paste0(ssfile, ".nii")
	############################
	# Copying files over
	############################
}

