rm(list=ls())
library(cttools)
library(fslr)
library(plyr)
library(extrantsr)
library(drammsr)
homedir = "~/"
rootdir = "~/CT_Registration/"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
rootdir = path.expand(rootdir)
homedir = file.path(rootdir, "Rorden_data")
progdir = file.path(rootdir, "programs")
basedir = file.path(rootdir, "Registration")

# scandir = file.path(homedir, "ct_scans")
regdir = file.path(homedir, "registered_ct_scans")
scandir = file.path(homedir, "ct_scans")
# setwd(homedir)
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")
outdir = file.path(homedir, "DRAMMS_Registered_Scans")

template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)
dtemp = dim(temp)

ss_template = file.path(tempdir, "scct_unsmooth_SS_0.01.nii.gz")
sstemp = readNIfTI(template)

ss_files = list.files(path=scandir, full.names=TRUE, 
	pattern="^.*_SS_0.01\\.nii")

outstubs = file.path(outdir, 
	paste0(nii.stub(ss_files, bn=TRUE), "_"))
reg_files = paste0(outstubs, "DRAMMS.nii.gz")

ifile = 1
for (ifile in seq_along(ss_files)){
	if (!file.exists(reg_files[ifile])){
		run = dramms(
			source = ss_files[ifile], 
			target = ss_template, 
			outfile = reg_files[ifile]
			)
	}
	print(ifile)
}

# invtransforms
# [1] "fileeb246cee66f30GenericAffine.mat" 
# [2] "fileeb246cee66f31InverseWarp.nii.gz"

