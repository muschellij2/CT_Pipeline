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


scandir = file.path(homedir, "CT_lesions")
# setwd(homedir)
outdir = tempdir()


files = list.files(scandir, pattern = "^v.*_CT.nii", full.names=TRUE)
ifile = 1

for (ifile in seq_along(files)){
	print(ifile)
	myfile = files[ifile]
	ssfile = paste0(nii.stub(myfile), "_SS_0.01")
	Sys.setenv(FSLOUTPUTTYPE = "NIFTI")
	ss = CT_Skull_Strip(myfile, 
		retimg=TRUE, 
		outfile = ssfile,
	 opts = "-f 0.01 -v")
	ssfile = paste0(ssfile, ".nii")
	############################
	# Copying files over
	############################
}

