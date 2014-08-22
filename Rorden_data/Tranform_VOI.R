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

scandir = file.path(homedir, "CT_lesions")
# setwd(homedir)
outdir = tempdir()


files = list.files(scandir, pattern = "cvp.*_sn.mat", full.names=TRUE)
ifile = 1


mid.folder = function(x, folname = ""){
  d = dirname(x)
  b = basename(x)
  file.path(d, folname, b)
}


for (ifile in seq_along(files)){
	print(ifile)
	mat = files[ifile]
	runfile = gsub("_CT_sn.mat$", "_voi.nii", mat)
	ssfile = file.path(  gsub("^c", "", basename(runfile)))
	runfile = file.path(dirname(runfile), 
		gsub("^c", "", basename(runfile)))
	ssmask = file.path(dirname(runfile),
		gsub("_voi", "_CT_SS_0.01_Mask", ssfile))
	ssfile = file.path(dirname(runfile),
		gsub("_voi", "_CT_SS_0.01", ssfile))


	run_ctnorm_write(mat, file=c(runfile, ssfile, ssmask))


}