##################################################################
## This code is for prediciton of Image Segmentation of CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
##################################################################
##################################################################
rm(list=ls())
library(plyr)
library(cttools)
library(ROCR)
library(matrixStats)
library(reshape2)
library(ggplot2)
library(fslr)
# library(car)
library(GGally)
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

segdir = file.path(progdir, "Segmentation")
source(file.path(segdir, "performance_functions.R"))

outdir = file.path(basedir, "results")
correct = "Rigid_sinc"
# options = c("none", "N3", "N4", "N3_SS", "N4_SS",
#         "SyN", "SyN_sinc", "Rigid", "Affine", "Rigid_sinc", 
#         "Affine_sinc")
options = c("none", "N3_SS", "N4_SS", 
		"Rigid", "Rigid_sinc")

# types = c("", "_include", "_zval", "_zval2")
types = c("_zval2", "_zval_all", '_zval2_medztemp')
# types = "_zval2"
# "_include_all", 
# types = "_include_all"
type = types[1]


eg = expand.grid(correct = options, 
	type = types, 
	stringsAsFactors=FALSE)
best.mats = best.worsts = NULL
vox.mats = vol.mats = NULL
group = "Test"

for (correct in options){

	print(correct)
	correct = match.arg(correct, options)
	adder = switch(correct, 
		"none"= "",
		"N3"="_N3",
		"N4" = "_N4",
		"N3_SS" = "_N3_SS",
		"N4_SS" = "_N4_SS", 
		"SyN" = "_SyN",
		"SyN_sinc" = "_SyN_sinc",
		"Rigid" = "_Rigid",
		"Affine" = "_Affine",
		"Rigid_sinc" = "_Rigid_sinc",
		"Affine_sinc" = "_Affine_sinc")

	for (type in types){
		
		outfile = file.path(outdir, 
			paste0("Best_", group, "_model_performance_results", 
				adder, type, 
				".Rda")
			)
		x = load(file=outfile)
		best.mat$group = group
		best.mat$correct = correct
		best.mat$type = gsub("^_", "" ,type)
		best.mats = rbind(best.mats, best.mat)

		best.worst$group = group
		best.worst$correct = correct
		best.worst$type = gsub("^_", "" ,type)
		best.worsts = rbind(best.worsts, best.worst)

		vox.eg$group = group
		vox.eg$correct = correct
		vox.eg$type = gsub("^_", "" ,type)

		vox.mats = rbind(vox.mats, vox.eg)

		vol.eg$group = group
		vol.eg$correct = correct
		vol.eg$type = gsub("^_", "" ,type)

		vol.mats = rbind(vol.mats, vol.eg)
		rm(list=x)
		print(type)
	} # type loop
} # correct loop

all.mats = rbind(vox.mats, vol.mats)

best.mats$volume = grepl("diff", best.mats$measure)
best.worsts$volume = grepl("diff", best.worsts$measure)

best.mats$sign = sign(best.mats$best.value)
best.worsts$sign = sign(best.worsts$best.worst)

best.mats$best.value[ best.mats$volume ] = 
	- abs(best.mats$best.value[ best.mats$volume ])

best.worsts$best.worst[ best.worsts$volume ] = 
	- abs(best.worsts$best.worst[ best.worsts$volume ])

best.mat = ddply(best.mats, .(measure), function(x){
	best = which.max(x$best.value)
	x$best.value = x$best.value * x$sign
	x[best, c("obj", "best", "best.value", "correct", 
		"type", "volume", "sign")]
})

best.worst = ddply(best.worsts, .(measure), function(x){
	best = which.max(x$best.worst)
	x$best.worst = x$best.worst * x$sign
	x[best, c("obj", "best", "best.worst", "correct", 
		"type", "volume", "sign")]
})

best.mat$best.value = 
	abs(best.mat$best.value) * best.mat$sign

best.worst$best.worst = 
	abs(best.worst$best.worst) * best.worst$sign 	

best.mat = best.mat[order(best.mat$volume, 
	best.mat$measure), ]
best.worst = best.worst[order(best.worst$volume, 
	best.worst$measure), ]


best.mats$best.value = 
	abs(best.mats$best.value) * best.mats$sign

best.worsts$best.worst = 
	abs(best.worsts$best.worst) * best.worsts$sign

best.mat$volume = best.worst$volume = NULL
best.mat$sign = best.worst$sign = NULL
best.mats$sign = best.worsts$sign = NULL

best.mat$mod = paste(best.mat$best, best.mat$correct, 
	best.mat$type, 
	sep= "_")

table(best.mat$mod)