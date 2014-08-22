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
regdir = file.path(homedir, "registered_ct_scans")


scandir = file.path(homedir, "ct_scans")
# setwd(homedir)
outdir = tempdir()


files = list.files(scandir, pattern = "*.nii.gz", full.names=TRUE)
ifile = 1

for (ifile in seq_along(files)){
	print(ifile)
	myfile = files[ifile]
	runfile = file.path(outdir, basename(myfile))
	file.copy(myfile, runfile, overwrite=TRUE)
	gunzip(runfile, overwrite = TRUE)
	runfile = gsub("\\.gz", "", runfile)

	code = "addpath(genpath('/home/bst/student/jmuschel/spm8')); "
	code = paste0(code, sprintf("clinical_c2h('%s');", runfile))
	run_matlab_code(code)
	outfile = file.path(dirname(runfile), 
		paste0("h", basename(runfile)))

	run_ctnorm(rawfile = outfile, deleteinter=FALSE)

	ssfile = paste0(nii.stub(outfile), "_SS_0.01")
	ss = CT_Skull_Strip(outfile, 
		retimg=TRUE, 
		outfile = ssfile,
	 opts = "-f 0.01 -v")
	gunzip(paste0(ssfile, ".nii.gz"))
	ssfile = paste0(ssfile, ".nii")
	############################
	# Copying files over
	############################
	file.gzcopy = function(from, todir, ...){
		file.copy(from, todir)
		g = file.path(todir, basename(from))
		gzip(g, overwrite=TRUE)
	}
	stub = nii.stub(basename(outfile))
	warped = file.path(outdir, 
		paste0("w", stub, ".nii"))
	file.gzcopy(warped, regdir, overwrite = TRUE)

	wmat = file.path(outdir, 
		paste0("c", stub, "_sn.mat"))
	file.copy(wmat, regdir, overwrite = TRUE)

	mat = file.path(outdir, 
		paste0(stub, ".mat"))
	file.copy(mat, regdir, overwrite = TRUE)

	warped.ss = run_ctnorm_write(wmat, 
		files= ssfile)

	ss.stub = nii.stub(basename(ssfile))
	sswarped = file.path(outdir, 
		paste0("w", ss.stub, ".nii"))
	file.gzcopy(sswarped, regdir, overwrite = TRUE)
}

