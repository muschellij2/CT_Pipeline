##########################################
## This code is for Creating Predictors 
# for each individual
##
## Author: John Muschelli
##########################################
##########################################
rm(list=ls())
library(ichseg)
library(methods)
library(neurobase)
library(dplyr)
rootdir = file.path("/Volumes/DATA_LOCAL", 
	"Image_Processing")
if (Sys.info()[["user"]] %in% "jmuschel") {
  rootdir = Sys.getenv("dex")
}
basedir = file.path(rootdir, 
	"PITCH_reconverted/processed_data")

ids = list.dirs(recursive = FALSE,
	path = basedir,
	full.names = TRUE)
imgs = sapply(ids, function(x) {
	f = list.files(pattern = "ROI.nii.gz", 
		path = x, 
		recursive = FALSE,
		full.names = TRUE)
})
df = data_frame(
	id = basename(ids),
	id_dir = ids,
	roi = imgs,
	img = sub("ROI[.]", ".", roi),
	stub = nii.stub(img, bn = TRUE))
df = df %>% 
	mutate(
		outdir = file.path(id_dir, 
			"processed"),
		outfile = file.path(outdir,
			"predictors.rds"))
iimg <- as.numeric(
	Sys.getenv("SGE_TASK_ID"))
if (is.na(iimg)) iimg = 111

runx = df[iimg,]
outdir = runx$outdir
outfile = runx$outfile

if (!file.exists(outfile)) {

	dir.create(outdir, recursive = TRUE,
		showWarnings = FALSE)
	img = runx$img
	roi = runx$roi
	stub = runx$stub
	if (!file.exists(roi)) {
		roi = NULL
	}
	proc = ich_process_predictors(
		img = img,
		n4_correct = FALSE,
		outdir = outdir,
		erode_mask = TRUE,
		roi = roi,
		save_imgs = TRUE,
		stub = stub)

	saveRDS(proc$img.pred, 
		outfile)  
}