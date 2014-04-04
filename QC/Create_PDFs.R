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


#### baseline NIHSSS ata 
nihss = read.csv(file.path(basedir, "baseline_NIHSS.csv"), 
                 stringsAsFactors=FALSE)
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")


outdir = file.path(basedir, "QC")
whichdir = "reoriented"

# whichdir = "FLIRT"
rerun = TRUE

template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)

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


view.pdf = function(fname){
	system(sprintf("xpdf %s&", fname))
}

view.png = function(fname){
	system(sprintf("display %s&", fname))
}

dtemp = dim(temp)
iid = 1;

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

	pdfname = file.path(outdir, paste0(id, ".pdf"))


	if (rerun | !file.exists(pdfname) ){

		# if (verbose) print(rawimgs[iimg])

		rawimg = readNIfTI(rawimgs[iimg])
		rawimg[is.nan(rawimg) | is.na(rawimg)] = 0
		# rawimg[ rawimg > 100 | rawimg < 0] = 0
		rawmask = rawimg >= 0 & rawimg <= 100

		# if (verbose) print(imgs[iimg])

		img = readNIfTI(imgs[iimg])
		img[is.nan(img) | is.na(img)] = 0
		tt = temp
		tt[ tt >= 99 | tt < 0] = 0
		tt[ img > 0 ] = 100

		## get cross hair information
		chair = round(colMeans(which(img > 0, arr.ind=TRUE)))

		chair[is.nan(chair)] = dtemp[is.nan(chair)]/2

		get.mar = function(img, margin=3){
			x = which(img > 0, arr.ind=TRUE)
			x = x[,margin]
			sort(unique(x))
		}
		zs = get.mar(img)

		rawimg = window_img(rawimg, replace="missing")
		img[img <= 0] = NA
		x = rawimg
		y = img

	    X <- nrow(x)
	    Y <- ncol(x)

		get.breaks = function(x, col.x){
			zlim.x = range(x, na.rm=TRUE)
	    	breaks.x <- c(min(x, zlim.x, na.rm = TRUE), 
	    		seq(min(zlim.x, na.rm = TRUE), max(zlim.x, na.rm = TRUE), 
	    			length = length(col.x) - 1), 
	    		max(x, zlim.x, na.rm = TRUE))
	    	return(list(zlim.x=zlim.x, breaks.x = breaks.x))
		}

		axes = FALSE
		plane = "axial"
		aspect <- x@pixdim[3]/x@pixdim[2]
		xlab = ylab = ""

		nzs = length(zs)
		pb = txtProgressBar(min=0, max=nzs, style=3)
		# pdf(pdfname, height= 3.5, width=7)
		pngs = file.path(outdir, paste0(id, seq(nzs), ".png"))

		for (iz in seq(nzs)){
		    z = zs[iz]

		    png(pngs[iz], 
		    	height= 3.5, 
		    	width=7, 
		    	units = "in",
		    	res= 300,
		    	type = "cairo")

			par(mfrow=c(1, 2), 
				oma = rep(0, 4), 
				mar = rep(0, 4), 
				bg = "black")
			### template image
			col.x = c(gray(0), gray(1:60/60), hotmetal(1))
			gb = get.breaks(tt, col.x)
			breaks.x = gb$breaks.x
			zlim.x = gb$zlim.x
            graphics::image(1:X, 1:Y, tt[, , z], 
            	col = col.x, 
                breaks = breaks.x, 
                zlim = zlim.x, 
                asp = aspect, 
                axes = axes, 
                xlab = xlab, 
                ylab = ylab)

            ### transformed patient information
			col.x = gray(0:64/64)
			gb = get.breaks(x, col.x)
			breaks.x = gb$breaks.x
			zlim.x = gb$zlim.x
			zlim.y <- range(y, na.rm = TRUE)		
			col.y = alpha("red", 0.25)

            graphics::image(1:X, 1:Y, x[, , z], col = col.x, 
                breaks = breaks.x, zlim = zlim.x, asp = aspect, 
                axes = axes, xlab = xlab, ylab = ylab)
            graphics::image(1:X, 1:Y, y[, , z], col = col.y, 
                zlim = zlim.y, add = TRUE)

            setTxtProgressBar(pb, value=iz)
			# overlay(rawimg, img, 
			# 	z=zs[1], plot.type = "single",
			# 	col.y=)
			dev.off()
        }
        close(pb)
		# dev.off()
		mystr= paste(pngs, collapse=" ", sep="")
		system(sprintf("convert %s %s", mystr, pdfname))
		file.remove(pngs)
	}

	print(iid)
}