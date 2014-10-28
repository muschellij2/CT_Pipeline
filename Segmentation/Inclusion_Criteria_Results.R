###################################################################
## This code is for prediction of Image Segmentation of CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
###################################################################
###################################################################
rm(list=ls())
library(plyr)
library(cttools)
library(fslr)
library(ROCR)
library(matrixStats)
library(reshape2)
library(ggplot2)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

outdir = file.path(basedir, "results")

correct = "N3_SS"
options = c("none", "N3", "N4", "N3_SS", "N4_SS",
        "SyN", "SyN_sinc", "Rigid", "Affine")
icorr <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(icorr)) icorr = 1
correct = options[icorr]

keep.obj = ls()


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
    "Affine" = "_Affine")


#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )

outfile = file.path(outdir, "111_Filenames.Rda")
load(file = outfile)


##############################
# Keeping files where predictors exist
##############################
outfiles = nii.stub(basename(fdf$img))
outfiles = paste0(outfiles, "_predictors", adder, ".Rda")
outfiles = file.path(fdf$outdir, outfiles)
stopifnot(file.exists(outfiles))

# load(file = file.path(outdir, "Segmentation_Models.Rda"))
##############################
# Run lmod number of models - not all the models - leave out
##############################
# mod.filename = file.path(outdir, 
# 	paste0("Collapsed_Models", adder, ".Rda"))
# load(mod.filename)

get.pred = 1
# get.pred = 50

fname = file.path(outdir, 
    paste0("Aggregate_data_cutoffs", adder, ".Rda"))

load(file = fname)
keepnames = colnames(est.cutoffs)
keepnames = keepnames[!( keepnames %in% 
	c("dist_centroid", "smooth10", "smooth20", "moment2", 
		"moment4", "moment3"))]

mod.filename = file.path(outdir, 
		paste0("Collapsed_Models", adder, ".Rda"))
load(mod.filename)	

nopred = seq(non.aggmods)
vd = vol.data
vd = vd[-nopred,]
nr = nrow(vd)
subset.ind = ceiling(nr/2)
subset.ind = seq(1, subset.ind)

ffdf = fdf[-nopred, ]
ffdf = ffdf[subset.ind,]

allnames = c(paste0(keepnames, ".cutoff"), 
		"include", "include.all", "zval", "zval2")
cn = c(outer(allnames, c(".out", ".reduced"), paste0))
ffdf[, cn] = NA
ffdf$nroi = NA

for (get.pred in seq(nrow(ffdf))){

	iddir = ffdf$iddir[get.pred]
	id.outdir = ffdf$outdir[get.pred]
	predname = nii.stub(basename(ffdf$img[get.pred]))
	predname = file.path(id.outdir, 
		paste0(predname, "_predictors", adder, ".Rda"))
	load(predname)
	df = img.pred$df
	stopifnot(all(df$Y %in% c(0, 1)))
	nim = img.pred$nim
	keep.ind = img.pred$keep.ind 
	rm(list="img.pred")
    for (i in 1:3) gc()	
	df$include = df$value >= 30 & df$value <= 100


    include = rep(TRUE, length=nrow(df))
    for (icut in keepnames){
    	qcuts = est.cutoffs[, icut]
    	colname = paste0(icut, ".cutoff")
    	df[, colname] = df[, icut] >= qcuts[1] & 
    		df[, icut] <= qcuts[2]
    	include = include & df[, colname]
    }

    df$zval = df[, "zscore3.cutoff"] & df$include &
    	df$pct_thresh.cutoff
    df$zval2 = df[, "zscore2.cutoff"] & df$zval

    # pdfname = file.path(outdir, 
    #     paste0("Aggregate_Data_Plots", adder, ".pdf"))
    # # pdf(pdfname)
    # pdfobj = smallpdf()
        
	df$include.all = include

	# qdat = melt(quants)
 #    colnames(qdat) = c("q", "pred", "value", "Y")
 #    qdat = qdat[ qdat$q %in% c("0.1%", "99.9%"), ]
 #    qdat$Y = as.numeric(qdat$Y)
 #    qdat = qdat[ qdat$Y == 1, ]
    # dd = df[ df$Y == 1, ]
    # h = ggplot(dd, aes(fill = factor(include.all))) + 
    #     geom_histogram()
    # iname = keepnames[1]

    # for (iname in keepnames){
    #     dat = qdat[ qdat$pred %in% iname, ]
    #     hh = h + aes_string(x = iname) + 
    #         geom_vline(data=dat, 
    #             aes(xintercept = value))
    #     hh = hh + ggtitle(iname)
    #     print(hh)
    #     print(iname)
    # } 



	df$in0100 = df$value >= 0 & df$value <= 100
	# df$in20_85 = df$value >= 20 & df$value <= 85
	df$mask = df$mask > 0
	
	######################################
	# Get volume of ROI
	######################################	
	not0100 = sum(df$Y[ !df$in0100 ])

	######################################
	# Keep all ROI = 1, even if not inmask
	######################################	
	roi.not.in = which(df$Y == 1)
	roi.not.in = roi.not.in[!(roi.not.in %in% keep.ind)]
	keep.ind = sort(c(keep.ind, roi.not.in))


	#### need this because the length of df has changed
	roi.not.in = which(keep.ind %in% roi.not.in)
	nroi.not.in = length(roi.not.in)

	df = df[keep.ind,]

	nroi = sum(df$Y == 1)

	pcts = matrix(NA, nrow=length(allnames), ncol=2)
	colnames(pcts) = c("out", "reduced")
	rownames(pcts) = allnames
	icut = allnames[1]
	for (icut in allnames){
		include = df[, icut]
		pct.out = sum(df$Y[! include]) / sum(df$Y)
		pct.reduced = sum(1-df$Y[!include]) / sum(1-df$Y)
		pcts[icut,] = c(pct.out, pct.reduced)
    }

    print(pcts)
    print(get.pred)
	predname = nii.stub(basename(ffdf$img[get.pred]))
	predname = file.path(id.outdir, 
		paste0(predname, "_cutoff_results", adder, ".Rda"))	
	save(pcts, allnames, keepnames, nroi,
		file = predname)

	df = melt(pcts)
	vals = df$value
	cn = paste0(df$Var1, ".", df$Var2)
	ffdf[get.pred, cn] = vals
	ffdf$nroi[get.pred] = nroi
} 

########################################
# Collapse results
########################################
# ffdf$nroi = NA
# for (get.pred in seq(nrow(ffdf))){

# 	iddir = ffdf$iddir[get.pred]
# 	id.outdir = ffdf$outdir[get.pred]
	
# 	predname = nii.stub(basename(ffdf$img[get.pred]))
# 	predname = file.path(id.outdir, 
# 		paste0(predname, "_cutoff_results", adder, ".Rda"))	
# 	load(file = predname)
# 	# ffdf$pct.out[ get.pred ] = pct.out
# 	# ffdf$pct.30out[ get.pred ] = pct.30out
# 	# ffdf$pct.reduced[ get.pred ] = pct.reduced
# 	# ffdf$pct.30reduced[ get.pred ] = pct.30reduced
# 	ffdf$nroi[get.pred] = nroi
# } 

predname = file.path(outdir, 
	paste0("Aggregate_cutoff_results", adder, ".Rda"))
cmeans = colMeans(ffdf[, cn])
cmaxs = colMaxs(ffdf[, cn])
cmins = colMins(ffdf[, cn])

x = data.frame(mean=t(t(cmeans)), max=t(t(cmaxs)), 
	min = t(t(cmins)))
x$var=  rownames(x)
x$cut = gsub("(.*)[.](out|reduced)", "\\1", x$var)
x$type = gsub("(.*)[.](out|reduced)", "\\2", x$var)
x$cut = gsub(".cutoff", "", x$cut)
rownames(x) = NULL
x$var = NULL
means = reshape(x, direction = "wide", idvar = c("cut"), 
	timevar = "type")

save(ffdf, keepnames, 
	cn, means,
	file = predname)

