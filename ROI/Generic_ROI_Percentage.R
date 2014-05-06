#################################
# Creates map of top p-values from unadusted model
# Author: John Muschelli
#################################
rm(list=ls())
library(oro.nifti)
library(plyr)
library(scales)
library(RColorBrewer)
library(data.table)
library(cttools)
library(fslr)
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
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

outdir = file.path(basedir, "results")

whichdir = "reoriented"
outcome = "NIHSS"
adder = paste0(outcome, "_")
if (outcome == "NIHSS"){
	adder = ""
}


rerun = FALSE

demog = read.csv(file=file.path(basedir, "Demog_NIHSS_Mask.csv"), 
	stringsAsFactors=FALSE)
demog$Base_ICH_10 = demog$Diagnostic_ICH /10

#### load voxel data
outfile = file.path(outdir, "Voxel_Matrix.Rda")
vmat = load(file=outfile )



load(file.path(atlasdir, "All_FSL_Atlas_Labels.Rda"))


usedf = jhut1.df
uselist = jhut1.list

get.ind = function(df, pattern, lrsub = TRUE){
  df$lab = tolower(df$Label)
  if (lrsub) {
    df$side = gsub("(.*)_(left|right)$", "\\2", df$lab)
    df$lab = gsub("(.*)_(left|right)$", "\\1", df$lab)
  }
  inds = df[grep(pattern, df$lab),]
  inds
}

####################################################################
# Find some pattern you want - like internal capsule
####################################################################
mypat = "internal_capsule"
inds = get.ind(usedf, pattern=mypat)

####################################################################
# Get the list of indices 
####################################################################
img.ind = uselist[inds$Label]

####################################################################
# Get number of voxels engaged for each pt
####################################################################
perc_engage = function(ind.list){
	sums = sapply(ind.list, function(x){
		submat = mat[x,]
		cs = colSums(submat)
	})

	####################################################################
	# Get total N of voxels for denominators
	lengths = sapply(ind.list, length)
	denom = matrix(lengths, nrow=nrow(sums), 
		ncol=length(ind.list),
		byrow=TRUE)
	#Percentage
	perc = sums/ denom
	return(list(sums=sums, perc = perc))
}

eng.sub.lr = perc_engage(img.ind)

####################################################################

####################################################################
# Collapse data by common label - shedding left/right designation

cimg.ind = dlply(inds, .(lab), function(x) {
	xx = x$Label
	u = unlist(img.ind[xx])
	u
})
eng.sub = perc_engage(cimg.ind)
####################################################################

####################################################################
# Collapse data by left/right
simg.ind = dlply(inds, .(side), function(x) {
	xx = x$Label
	u = unlist(img.ind[xx])
	u
})
eng.lr = perc_engage(simg.ind)


xx = inds$Label
all.ind = list(unlist(img.ind[xx]))
names(all.ind) = mypat

eng.all = perc_engage(all.ind)


cor(demog$Q6a_motor_left_leg, eng.sub$perc, use="na.or.complete")
cor(demog$Q6a_motor_left_leg, eng.lr$perc, use="na.or.complete")
cor(demog$Q6a_motor_left_leg, eng.sub.lr$perc, use="na.or.complete")
cor(demog$Q6a_motor_left_leg, eng.all$perc, use="na.or.complete")


cor(demog$Q6b_motor_right_leg, eng.sub$perc, use="na.or.complete")
cor(demog$Q6b_motor_right_leg, eng.lr$perc, use="na.or.complete")
cor(demog$Q6b_motor_right_leg, eng.sub.lr$perc, use="na.or.complete")
cor(demog$Q6b_motor_right_leg, eng.all$perc, use="na.or.complete")

