```{r knit-setup, echo=FALSE, results='hide'}
library(knitr)
opts_chunk$set(message=FALSE, warning=FALSE,
               comment="",
               echo=FALSE, results='hide')
opts_knit$set(error = FALSE)

# options(rstudio.markdownToHTML =
#   function(inputFile, outputFile) {
#     library(knitrBootstrap)
#     knit_bootstrap_md(input=inputFile, output=outputFile, code_style="Brown Paper", chooser=c("boot", "code"), show_code=FALSE)
#   }
# )
```
 

```{r, echo=FALSE}
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
### need cairo for cluster
options(bitmapType='cairo')
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
source(file.path(progdir, "convert_DICOM.R"))
source(file.path(progdir, "fslhd.R"))
basedir = file.path(rootdir, "Registration")


#### baseline NIHSSS ata 
nihss = read.csv(file.path(basedir, "baseline_NIHSS.csv"), 
                 stringsAsFactors=FALSE)
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

outdir = file.path(basedir, "results")
whichdir = "reoriented"

# whichdir = "FLIRT"
rerun = FALSE

setwd(progdir)


all.ids = list.dirs(basedir, recursive=FALSE, full.names=FALSE)
all.ids = all.ids[grepl("\\d\\d\\d-(\\d|)\\d\\d\\d", all.ids)]
redir = file.path(basedir, all.ids, whichdir)
redir = redir[file.exists(redir)]
nfiles = sapply(redir, function(x) 
	length(dir(path=x, pattern="ROI*.nii.*")))
redir = redir[nfiles > 0]
## those that have ROIs
all.ids = gsub(paste0(basedir, "/(.*)/", whichdir), "\\1", redir)
# all.ids = all.ids[file.exists(redir)]


```{r}

ID = "100-365"
run.id = function(ID, whichdir, rerun=TRUE, verbose=TRUE){
	### go to id directory
	# if (verbose) print(ID)
	iddir = file.path(basedir, ID)
	redir = file.path(iddir, whichdir)
	if (!file.exists(redir)){
		system(sprintf('mkdir -p "%s"', redir))
	}
	imgs = list.files(pattern='ROI.*.nii', full.names=TRUE, path=redir)
	imgs = imgs[ !grepl("2mm_", imgs)]
	imgs = imgs[ !grepl("affine9", imgs)]
	imgs = imgs[ !grepl("bws", imgs)]

	# rawimgs = gsub("bws", "w", imgs)
	rawimgs = gsub("^ROI_", "", imgs)
	rawimgs = gsub("ROI\\.", ".", rawimgs)

	# rois = lapply(imgs, readNIfTI)
	keep.imgs = 1
	nimgs = length(imgs)
	ids = pngs = rep(NA, nimgs)


	resdir = file.path(redir, "results")
	if (!file.exists(resdir)){
		system(sprintf('mkdir -p "%s"', resdir))
	}
	iimg = 1;


	# for (iimg in seq(nimgs)){
		# if (verbose) print(iimg)

		id = basename(rawimgs[iimg])
		id = gsub("\\.gz$", "", id)
		iid = id = gsub("\\.nii$", "", id)
		id = gsub("^w", "spm_", id)
		id = gsub("^affine", "flirt_dof", id)

		ids[iimg] = gsub("^(w|affine(9|12)_)", "", iid)
		pngs[iimg] = file.path(resdir, paste0("native_", id, ".png"))

		if (rerun | 
			!file.exists(pngs[iimg]) ){

			if (verbose) print(rawimgs[iimg])
			if (!file.exists(rawimgs[iimg]) | 
				!file.exists(imgs[iimg]) ) return(NULL);
			rawimg = readNIfTI(rawimgs[iimg], reorient=FALSE)
			rawimg[is.nan(rawimg) | is.na(rawimg)] = 0
			# rawimg[ rawimg > 100 | rawimg < 0] = 0
			rawmask = rawimg >= 0 & rawimg <= 100

			if (verbose) print(imgs[iimg])

			img = readNIfTI(imgs[iimg], reorient=FALSE)
			img[is.nan(img) | is.na(img)] = 0

			## get cross hair information
			chair = round(colMeans(which(img > 0, arr.ind=TRUE)))

			if (any(is.na(chair))) chair= NULL
			png(pngs[iimg], type="cairo")

				mask.overlay(rawimg, img, xyz=chair, 
					col.y=alpha("red", 0.25))

			dev.off()

		}
	# }


	return(list(nimgs=nimgs, 
		ids=ids,
		pngs=pngs))
}

# run1 = run.id("100-318", whichdir=whichdir, rerun=TRUE)
res = llply(all.ids, run.id, whichdir=whichdir, rerun=rerun,
	.progress ="text")
```
