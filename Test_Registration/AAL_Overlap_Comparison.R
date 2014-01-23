#####################################
## Author: John Muschelli
## Date: January 20, 2014
## Purpose: Read in the AAL atlas labels and make R objects
## that can be used later for overlap metrics.
#####################################
#####################################
rm(list=ls())
library(R.matlab)
library(oro.nifti)
library(plyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
homedir = "/Applications"
basedir = "/Volumes/DATA_LOCAL/Image_Processing/Test_Registration"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  basedir = "/dexter/disk2/smart/stroke_ct/ident/Test_Registration"
}
spmdir = file.path(homedir, "spm8")
aaldir = file.path(spmdir, "toolbox", "aal_for_SPM8")
load(file=file.path(aaldir, "ROIList.Rda"))

load(file=file.path(basedir, "FLIRT", 
	"ROI_Overlap_Measures.Rda"))
flirt = allres
flirt$col_name = NULL

load(file=file.path(basedir, "reoriented", 
	"ROI_Overlap_Measures.Rda"))
spm = allres
spm$col_name = NULL

flirt.ids = unique(flirt$fname)
spm.ids = unique(spm$fname)

sdiff = setdiff(flirt.ids, spm.ids)
stopifnot(length(sdiff) == 0)

sdiff = setdiff(spm.ids, flirt.ids)
stopifnot(length(sdiff) == 0)

allres = merge(spm, flirt, 
	by=c("fname", "area"), 
	suffixes = c(".spm", ".flirt"),
	sort=FALSE)
rois$area = rois$name
allres = merge(allres, rois[, c("area", "col_name")],
	by="area", all.x=TRUE,
	sort=FALSE)

nas = sapply(allres, is.na)
stopifnot(!any(nas))

#### drop any region where it was never categorized in either method
no.data = allres$raw.flirt == 0 & allres$raw.spm == 0

allres = allres[! no.data, ]

### make a long dataset for plotting/tabulation
long = melt(allres, id.vars=c("fname", "area", "col_name"))
## use collapsed name
long$area = NULL
long = aggregate(value ~ ., data=long, sum)
long$area = long$col_name
long$col_name = NULL
long$variable = as.character(long$variable)
ss = strsplit(long$variable, "\\.")
long$type = sapply(ss, function(x) x[1])
long$proc = sapply(ss, function(x) x[2])

ss = strsplit(long$fname, "_")
long$id = sapply(ss, function(x) x[1])

tmp = long[ long$type == "raw" & long$proc == "spm", ]

colourCount = length(unique(tmp$area))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

g = qplot(data=tmp, x=id, fill=factor(area), y=value,
	geom="bar", 
	stat="identity") +  
	scale_fill_manual(values = getPalette(colourCount)) +coord_flip()

# areas = tmp[ tmp$area %in% c("Caudate", "Putamen", "Thalamus")]
# g = qplot(data=tmp, x=id, fill=factor(area), y=value,
# 	geom="bar", 
# 	stat="identity") +  
# 	scale_fill_manual(values = getPalette(colourCount)) +coord_flip()



	# cres = cres[order(cres$weighted),]
# }