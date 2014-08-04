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



loadfile = file.path(atlasdir, "All_FSL_Atlas_Labels.Rda")
load(loadfile)

jhut1.col.df = jhut1.nolr.df
jhut1.col.df$Orig = jhut1.col.df$Label
jhut1.col.df$Label = tolower(jhut1.col.df$Label)
jhut1.col.df$Label = gsub("^(posterior|superior|inferior|anterior)_", 
	"", jhut1.col.df$Label)
	jhut1.col.df$Label = gsub("^(middle|lateral)_", 
	"", jhut1.col.df$Label)
jhut1.col.df$Location = NA
jhut1.col.df$Location[ grepl( "thalam",  jhut1.col.df$Label)] = 
	"Thalamus"
jhut1.col.df$Location[ grepl( "puta",  jhut1.col.df$Label)] = 
	"Putamen"
jhut1.col.df$Location[ grepl( "obus",  jhut1.col.df$Label)] = 
	"Globus Pallidus"

jhut1.col.df$Location[ grepl( "occipital|temporal",  jhut1.col.df$Label)] = 
	"Lobar"	
jhut1.col.df$Location[ grepl( "occipital",  jhut1.col.df$Label)] = 
	"Lobar"		
jhut1.col.df$Location[ grepl( "frontal_gyrus",  jhut1.col.df$Label)] = 
	"Lobar"		
jhut1.col.df$Location[ grepl( "fronto-orbital",  jhut1.col.df$Label)] = 
	"Lobar"		
jhut1.col.df$Location[ grepl( "frontal_wm",  jhut1.col.df$Label)] = 
	"Lobar"			
jhut1.col.df$Location[ grepl( "parietal_wm",  jhut1.col.df$Label)] = 
	"Lobar"			

jhut1.col.df$Location[ grepl( "cerebel",  jhut1.col.df$Label)] = 
	"Cerebellum"	

jhut1.col.df$Location[ grepl( "midbrain|pons|medulla",  jhut1.col.df$Label)] = 
	"Brain Stem"		

jhut1.col.df[ is.na(jhut1.col.df$Location),]




no.uncat = TRUE
for (no.uncat in c(TRUE, FALSE)){
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

	dd2 = merge(dd, max.area, all=TRUE, by="id")
	dd2$cat.area = dd2$area
	dd2$cat.area[grep("Lobe", dd2$cat.area)] = "Lobar"


	print(with(dd2, table(cat.area, clot_location)))
}