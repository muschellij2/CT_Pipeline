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
library(smallpdf)
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


#### baseline NIHSSS ata 
nihss = read.csv(file.path(basedir, "baseline_NIHSS.csv"), 
                 stringsAsFactors=FALSE)
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")


outdir = file.path(basedir, "BrainQC")
whichdir = "reoriented"

# whichdir = "FLIRT"
rerun = TRUE

template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)
tmask = file.path(tempdir, "scct_mask.nii.gz")
tempmask = readNIfTI(tmask)

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
uid = as.numeric(gsub("-", "", all.ids))

### get clot location
demog = read.csv(file.path(basedir, 
	"All_180_FollowUp_wDemographics.csv"), 
                 stringsAsFactors=FALSE)
demog$Sex = factor(demog$Gender, levels=c("Female", "Male"))
demog$Diagnostic_ICH = demog$ICH_Dx_10 *10

demog = demog[ demog$patientName %in% uid, ]
dd = demog[ order(demog$patientName), ]
dd = dd[, c("patientName", "Clot_Location_RC")]
colnames(dd) = c("id", "clot_location")


dtemp = dim(temp)
iid = 1;

tt = tempmask
tt[ tt < 1] = NA

## get cross hair information

get.mar = function(img, margin=3){
	x = which(img > 0, arr.ind=TRUE)
	x = x[,margin]
	sort(unique(x))
}
zs = get.mar(tempmask)
nzs = length(zs)



for (iid in seq_along(all.ids)){
	

	ID = all.ids[iid]
	iddir = file.path(basedir, ID)
	redir = file.path(iddir, whichdir)

	imgs = list.files(pattern='ROI.*.nii', full.names=TRUE, path=redir)
	imgs = imgs[ !grepl("2mm_", imgs)]
	if (grepl("reoriented", whichdir)) imgs = imgs[ grepl("bws", imgs)]
	imgs = imgs[ !grepl("affine9", imgs)]

	rawimgs = gsub("bws", "w", imgs)
	rawimgs = gsub("^ROI_", "", rawimgs)
	rawimgs = gsub("ROI\\.", ".", rawimgs)

	bin.popimg = popimg = array(0, dim=dtemp)

	# rois = lapply(imgs, readNIfTI)
	keep.imgs = 1
	imgs = imgs[keep.imgs]
	rawimgs = rawimgs[keep.imgs]
	nimgs = length(imgs)
	ids = opngs = rawpngs = rep(NA, nimgs)

	iimg = 1;
	# for (iimg in seq(nimgs)){
		# if (verbose) print(iimg)

	id = basename(rawimgs[iimg])
	id = nii.stub(id)
	id = gsub("^w", "", id)

	pdfname = file.path(outdir, paste0(id, "_Template.pdf"))

	if (rerun | !file.exists(pdfname) ){

		rawimg = readNIfTI(rawimgs[iimg])
		rawimg[is.nan(rawimg) | is.na(rawimg)] = 0
		rawimg[ rawimg > 100 | rawimg < 0] = 0
		rawimg = window_img(rawimg)
		# if (verbose) print(imgs[iimg])


		pb = txtProgressBar(min=0, max=nzs, style=3)

		sp = smallpdf()
		par(mfrow=c(1, 1), bg = "black")
		for (iz in seq(nzs)){
		    z = zs[iz]
			### template image

			overlay(rawimg, tt, z = z, plot.type = "single", 
				col.y=alpha("red", .5))
            setTxtProgressBar(pb, value=iz)
        }
        close(pb)
        smallpdf.off(pdfname, mypattern =sp$mypattern, dev = sp$dev)
	}

	print(iid)
}