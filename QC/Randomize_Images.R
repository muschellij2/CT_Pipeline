rm(list=ls())
library(R.matlab)
library(oro.nifti)
library(plyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(knitrBootstrap)
library(xtable)
library(scales)
library(fslr)
library(cttools)
### need cairo for cluster
options(bitmapType='cairo')
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
# source(file.path(progdir, "convert_DICOM.R"))
# source(file.path(progdir, "fslhd.R"))
basedir = file.path(rootdir, "Registration")


outdir = file.path(basedir, "QC")

set.seed(20140401)
pdfs = list.files(path=outdir, 
	pattern="\\d\\d\\d-\\d\\d\\d.*\\.pdf")
pdfs = data.frame(pdf =pdfs, stringsAsFactors=FALSE)
pdfs$id = gsub("^(\\d\\d\\d-\\d\\d\\d(|\\d))_.*", 
	"\\1", pdfs$pdf)

pdfs$nii = gsub("\\.pdf$", ".nii.gz", pdfs$pdf)
N = nrow(pdfs)
### 3 readers
readers = c("Natalie", "Sam", "Andrew")
fnames = file.path(outdir, pdfs$pdf)

nii.fnames = gsub("\\.pdf$", ".nii.gz", fnames)

### re-runs - duplicates - so we can get intra-reader consistency
perc = 0.2
n = round(N*perc)
dup.samp = pdfs[sample(n),]

for (iread in readers){
	pdfs[, iread] = sample(N, replace=FALSE)
	read.dir = file.path(outdir, iread)
	if (!file.exists(read.dir)){
		dir.create(read.dir)
	}
	new.fnames = file.path(read.dir, 
		sprintf("%04.0f.pdf", pdfs[, iread]))
	file.copy(fnames, new.fnames, overwrite=TRUE)

##### copying over nifti files
	new.nii = file.path(read.dir, 
		sprintf("%04.0f.nii.gz", pdfs[, iread]))
	file.copy(nii.fnames, new.nii, overwrite=TRUE)


	dup.samp[, iread] = sample(n, replace=FALSE) + N

	new.fnames = file.path(read.dir, 
		sprintf("%04.0f.pdf", dup.samp[, iread]))
	file.copy(file.path(outdir, dup.samp$pdf), 
		new.fnames, 
		overwrite=TRUE)

	new.nii = file.path(read.dir, 
		sprintf("%04.0f.nii.gz", dup.samp[, iread]))
	file.copy(file.path(outdir, dup.samp$nii), 
		new.nii, 
		overwrite=TRUE)

}

mydate = date()
save(dup.samp, perc, readers, pdfs, mydate,
	file=file.path(outdir, "Reader_Data.Rda"))


