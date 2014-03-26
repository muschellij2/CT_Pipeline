#################################
# Performs voxel-wise models for age/gender/volume
#################################
rm(list=ls())
library(oro.nifti)
library(plyr)
library(limma)
library(microbenchmark)
library(abind)
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
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

outdir = file.path(basedir, "results")

whichdir = "reoriented"

outfile = file.path(atlasdir, 
  paste0(whichdir, "_All_Atlas_ROI_Overlap_Measures.Rda"))

load(file=outfile)

makeres = function(allres){
	allres$fname = gsub("^bws", "", allres$fname)
	allres$fname = gsub("ROI$", "", allres$fname)
	allres[, c("raw", "bin", "weighted")] = 
		round(allres[, c("raw", "bin", "weighted")] *100, 1)
	ss = strsplit(allres$fname, "_")
	allres$id = sapply(ss, function(x) x[1])
	return(allres)
}
mni.allres = makeres(mni.allres)

uid = unique(as.numeric(gsub("-", "", mni.allres$id)))

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

no.uncat = TRUE

max.area = ddply(mni.allres, .(fname), function(x){
	if (no.uncat) x = x[!(x$area %in% "Uncategorized"),]
	ind = which(x$raw == max(x$raw))
	stopifnot(length(ind) == 1)
	x[ind, c("id", "area")]
})

max.area = ddply(max.area, .(id), function(x){
	x[1,]
})

max.area$id = as.numeric(gsub("-", "", max.area$id))
max.area = max.area[, c("id", "area")]

dd = merge(dd, max.area, all=TRUE, by="id")
dd$cat.area = dd$area
dd$cat.area[grep("Lobe", dd$cat.area)] = "Lobar"




















with(dd, table(cat.area, clot_location))