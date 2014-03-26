#################################
# Performs voxel-wise models for age/gender/volume
#################################
rm(list=ls())
library(oro.nifti)
library(plyr)
library(limma)
library(microbenchmark)
library(scales)
library(RColorBrewer)
library(data.table)
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

#### load voxel data
outfile = file.path(outdir, "Voxel_Matrix.Rda")
load(file=outfile )

outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile)

outfile = file.path(outdir, "Pvalue_Matrix.Rda")
load(file=outfile)


group=setkey(group, 'group')
ugroup = unique(group)
urows = ugroup$id
ugroups = ugroup$group
#### keeping if over 10 people have ICH in that locaiton



template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)
dtemp = dim(temp)

stopifnot(prod(dtemp) == nrow(mat))

thresholds = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1)
imod = 1;

recenter = TRUE

open.dev = function(file, type= "cairo", ...){
	get.ext = gsub("(.*)\\.(.*)$", "\\2", file)
	stopifnot( get.ext %in% c("pdf", "bmp", "svg", "png", 
		"jpg", "jpeg", "tiff"))

	## device is jpeg
	if (get.ext == "jpg") get.ext = "jpeg"
	### difff arguments for diff devices
	if (get.ext %in% c("pdf")) {
		do.call(get.ext, list(file=file, ...))
	} else if (get.ext %in% 
		c("bmp", "jpeg", "png", "tiff", "svg")) {
		do.call(get.ext, list(filename=file, type= type,...))
	}
}

device = "pdf"

view.png = function(fname){
	system(sprintf("display %s", fname))
}

view.pdf = function(fname){
	system(sprintf("xpdf %s", fname))
}

for (imod in seq(dim(res)[3])){


	fp = file.path(outdir, 
		paste0("Regression_Map_heatcol", imod, ".", device))
	open.dev(fp)
		res.p = temp
		res.p[!is.na(res.p)] = NA
		y = res[,"X",imod]
		cuts = cut(y, breaks=thresholds, include.lowest=TRUE)
		cuts = factor(cuts, levels=rev(levels(cuts)))
		cuts = as.numeric(cuts)-1
		# res.p[ rs > ncut ] = cuts
		res.p[rs > ncut] = -log(y)
		# res.p[rs > ncut] = as.numeric(y <= 0.05)
		# res.p[ rs > ncut ]  = 1-y
		cols = c("blue", "green", "yellow", "orange", "red")
		cols = alpha(cols, .7)
		mask.overlay(temp, res.p, col.y=cols)
		col.y = alpha(hotmetal(10), 0.7)	
		mask.overlay(temp, res.p, col.y=col.y)
		col.y[1] = alpha("white", 0.7)
		mask.overlay(temp, res.p, col.y=col.y)
		cols = alpha(brewer.pal(9, "Reds"), 0.7)
		mask.overlay(temp, res.p, col.y=cols)
		rm(list="res.p")
		
	dev.off()	


	for (recenter in c(TRUE, FALSE)){
		app = ""
		if (recenter) app = "_centered"


		### 
		fp = file.path(outdir, 
			paste0("Regression_Map", imod, 
				app, ".", device))
		open.dev(fp)
			res.p = temp
			res.p[!is.na(res.p)] = NA
			y = res[,"X",imod]
			cuts = cut(y, breaks=thresholds, include.lowest=TRUE)
			cuts = factor(cuts, levels=rev(levels(cuts)))
			cuts = as.numeric(cuts)-1
			# res.p[ rs > ncut ] = cuts
			res.p[rs > ncut] = as.numeric(y <= 0.05)
			xyz= NULL
			if (any(y <= 0.05) & recenter) {
				xyz = floor(colMeans(which(res.p > 0, arr.ind = TRUE)))
			}			
			# res.p[ rs > ncut ]  = 1-y
			# cols = c("blue", "green", "yellow", "orange", "red", "white")
			mask.overlay(temp, res.p, col.y="red", xyz= xyz)
			rm(list="res.p")			
		dev.off()
	# system(sprintf("xpdf %s", fp))


	#### bonferroni adjustment
		fp = file.path(outdir, 
			paste0("Regression_Map_Bonf_", imod, 
				app, ".", device))
		open.dev(fp)
			res.p = temp
			res.p[!is.na(res.p)] = NA
			y = res[,"X",imod]		
			# y = p.adjust(y, method="bonferroni", n=nuniq.rows)
			y = pmin(y * nuniq.rows, 1)
			if (any(y <= 0.05)) cat("Significant p-values\n")
			print(min(y))
			print(paste0("Min of Bonf corrected p-values: ", min(y)))
			res.p[rs > ncut] = as.numeric(y <= 0.05)		
			xyz= NULL
			if (any(y <= 0.05) & recenter) {
				xyz = floor(colMeans(which(res.p > 0, arr.ind = TRUE)))
			}		
			# res.p[ rs > ncut ]  = 1-y
			mask.overlay(temp, res.p, col.y="red", xyz= xyz)
			rm(list="res.p")

		dev.off()

	#### fdr adjustment
		fp = file.path(outdir, 
			paste0("Regression_Map_FDR_", imod, 
				app, ".", device))
		open.dev(fp)
			res.p = temp
			res.p[!is.na(res.p)] = NA
			y = res[,"X",imod]		
			y = p.adjust(y, method="fdr")
			if (any(y <= 0.05)) cat("Significant p-values\n")
			print(paste0("Min of FDR corrected p-values: ", min(y)))
			res.p[rs > ncut] = as.numeric(y <= 0.05)
			xyz= NULL
			if (any(y <= 0.05) & recenter) {
				xyz = floor(colMeans(which(res.p > 0, arr.ind = TRUE)))
			}		
			# res.p[ rs > ncut ]  = 1-y
			mask.overlay(temp, res.p, col.y="red", xyz= xyz)
			rm(list="res.p")
		dev.off()	

	#### fdr with diff p
		fp = file.path(outdir, paste0("Regression_Map_FDR_red_", 
			imod, app, ".", device))
		open.dev(fp)
			res.p = temp
			res.p[!is.na(res.p)] = NA
			y = res[,"X",imod]		
			y = y[urows]
			y = p.adjust(y, method="fdr")
			if (any(y <= 0.05)) cat("Significant p-values\n")
			print(paste0("Min of reduced FDR corrected p-values: ", 
				min(y)))
			y = data.frame(pval=y)
			y$group = ugroups
			dt.y = data.table(y, key='group')
			setkey(group, 'group')
			y = dt.y[group]
			stopifnot(!any(is.na(y$pval)))
			setkey(y, key='id')
			y = y$pval
			res.p[rs > ncut] = as.numeric(y <= 0.05)
			xyz= NULL
			if (any(y <= 0.05) & recenter) {
				xyz = floor(colMeans(which(res.p > 0, arr.ind = TRUE)))
			}

			# res.p[ rs > ncut ]  = 1-y
			mask.overlay(temp, res.p, col.y="red", xyz= xyz)
			rm(list="res.p")			
		dev.off()	

		print(imod)
	}

}