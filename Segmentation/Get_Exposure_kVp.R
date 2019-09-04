###################################################################
## This code is for Image Segmentation of CT
## The code is R but calls out FSL
##
## Author: John Muschelli
## Last updated: May 20, 2014
###################################################################
###################################################################
rm(list=ls())
library(plyr)
library(fslr)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
segdir = file.path(progdir, "Segmentation")

basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")
outdir = file.path(basedir, "results")

#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )

outfile = file.path(outdir, "111_Filenames_with_volumes_stats.Rda")
xxx = load(file = outfile)

iimg <- suppressWarnings(as.numeric(Sys.getenv("SGE_TASK_ID")))
if (is.na(iimg)) iimg = 34
## 2,12,17,34,37, 46, 48, 68, 70,81, 85, 86, 87, 99
#  has variable slice thickness
## 15 is not correct for sthickness
## 17 & 87 worst - has overlapping slice somewhat
## 75 has all negative
## 71 has no position data
## 13,71,101 has spacing

fdf$median = fdf$mode = fdf$mean = NA

fdf$thickvol = fdf$zvol = fdf$varslice = 
fdf$gantry = fdf$truevol = NA

# fdf = fdf[c(2,12,17,34,37, 46, 48, 68, 70,81, 85, 86, 87, 99),]
# dcmtables[, '0018-1152-Exposure']
alltabs = exposures = vector(mode="list", length=nrow(fdf))
vdims = exposures


for (iimg in seq(nrow(fdf))){
	
	runx = x = fdf[iimg,]
	sortdir = file.path(x$iddir, "Sorted")

	# run_model = function(x, fpr.stop = .1){
	fname = xfname = nii.stub(x$img, bn=TRUE)
	rda_stub = file.path(sortdir, fname)
	rda = paste0(rda_stub, "_ungantry_Header_Info.Rda")
	if (file.exists(rda)){
		# stop("OK here's one")
	} else {
		rda = paste0(rda_stub, "_Header_Info.Rda")
	}
	xrda = load(rda)

	alltabs[[iimg]] = dcmtables

	vdims[[iimg]] = voxdim(x$ssimg)
	# print(grep("pac", colnames(dcmtables), value=TRUE))
	# print(iimg)

	# cn = c("0018-5100-PatientPosition", 
	# 	"0020-0032-ImagePositionPatient")
	###############################
	# Slice thickness or Patient position?
	###############################
	gant = unique(as.numeric(
		dcmtables[, "0018-1120-GantryDetectorTilt"]))
	stopifnot(length(gant) == 1)

	cn = c("0020-0032-ImagePositionPatient")	
	dcmnames = colnames(dcmtables)
	expo = dcmtables[, 
		grepl("xposure", dcmnames) | 
			grepl("KVP", toupper(dcmnames)) |
			grepl("ROW", toupper(dcmnames)) |
			grepl("COLUMN", toupper(dcmnames)) |
			grepl("PixelSpacing", dcmnames) |
			grepl("SliceThickness", dcmnames)
			, 
		drop=FALSE]
	rownames(expo) = NULL
	exposures[[iimg]] = expo
}

ex = lapply(exposures, function(x){
	cn = colnames(x)
	last = function(r) r[length(r)]
	cn = sapply(strsplit(cn, "-"), last)
	colnames(x) = cn
	x
	})

vdims = do.call("rbind", vdims)
vsizes = apply(vdims, 1, prod)


st = as.numeric(
	unlist(sapply(sapply(ex, `[[`, "SliceThickness"), 
	unique)))
st = sort(unique(st))

ps = unique(
	unlist(sapply(sapply(ex, `[[`, "PixelSpacing"), 
	unique)))
ps = as.numeric(sapply(strsplit(ps, " "), `[`, 1))
ps = sort(unique(ps))


kvps = unique(unlist(sapply(sapply(ex, `[[`, "KVP"), unique)))
kvps = sort(kvps[!kvps %in% c("")])
kvps = kvps[!is.na(kvps)]
kvps = sort(as.numeric(kvps))

exps = sapply(ex, function(x){
	cn = colnames(x)
	if ("Exposure" %in% x){
			
		} else {
			return(NA)
		}
})

exps = unique(unlist(sapply(sapply(ex, `[[`, "Exposure"), unique)))
exps = exps[!exps %in% c("")]
exps = exps[!is.na(exps)]
exps = sort(as.numeric(exps))


scan.outfile = file.path(outdir, "Scanning_Parameters.Rda")
save(vdims, exps, vsizes, st, ps, kvps, ex,
	alltabs, fdf,
	file = scan.outfile)

