rm(list=ls())
library(cttools)
library(fslr)
library(plyr)
library(extrantsr)
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
lesdir = file.path(homedir, "CT_lesions")
# setwd(homedir)
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")
outdir = file.path(homedir, "SyN_Registered_Lesions")

template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)
dtemp = dim(temp)

ss_template = file.path(tempdir, "scct_unsmooth_SS_0.01.nii.gz")
sstemp = readNIfTI(template)

ss_files = list.files(path=lesdir, full.names=TRUE, 
	pattern="^.*_CT_SS_0.01\\.nii")

voi_files = gsub("CT_SS_0[.]01", "voi", ss_files)

outstubs = file.path(outdir, 
	paste0(nii.stub(ss_files, bn=TRUE), "_"))
reg_files = paste0(outstubs, "SyN.nii.gz")

voi.ofiles = file.path(outdir, 
	paste0(nii.stub(voi_files, bn=TRUE), "_SyN.nii.gz"))

ifile = 1
for (ifile in seq_along(ss_files)){
	outprefix = outstubs[ifile]
	if (!file.exists(reg_files[ifile]) |
		!file.exists(voi.ofiles[ifile])
		){
		run = ants_regwrite(
			filename = ss_files[ifile], 
			template.file = ss_template, 
			outfile = reg_files[ifile],
			other.files = voi_files[ifile],
			other.outfiles = voi.ofiles[ifile],
			outprefix = outprefix)
	}
	warp = paste0(outprefix, "1Warp.nii.gz")
	if (file.exists(warp)){
		file.remove(warp)
	}
	print(ifile)
}

# invtransforms
# [1] "fileeb246cee66f30GenericAffine.mat" 
# [2] "fileeb246cee66f31InverseWarp.nii.gz"

