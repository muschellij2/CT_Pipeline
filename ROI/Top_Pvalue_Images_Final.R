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
library(smallPDF)
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

#### load voxel data
outfile = file.path(outdir, "Voxel_Matrix.Rda")
vmat = load(file=outfile )

outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile)

###########################################
# Get 111 patients
###########################################

ids.111 = read.csv(file.path(basedir, "111_patients.csv"), 
	stringsAsFactors= FALSE)
uid = ids.111$patientName
all.ids = ids.111$id

mat.ids = colnames(mat)
mat.ids = gsub(".*(\\d\\d\\d-(\\d|)\\d\\d\\d).*", "\\1", mat.ids)
mat = mat[ , which(mat.ids %in% all.ids)]

#### reading in template
template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)
dtemp = dim(temp)


area_pct = function(img, ind.list, keepall) {
  ## get overlap of indices
  raw.mat = sapply(ind.list, function(x) sum(img[x]))
  any.mat = sapply(ind.list, function(x) mean(img[x] > 0))
  mn.mat = sapply(ind.list, function(x) mean(img[x]))
  names(raw.mat) = names(ind.list)
  ## cs is sum of indices of overlap
  cs.raw = data.frame(nvox=raw.mat, roi_pct_any = any.mat,
  	roi_mean_pct = mn.mat) 
  rownames(cs.raw) = names(ind.list)
  if (!keepall) cs.raw = cs.raw[rowSums(cs.raw) != 0, , drop=FALSE]
  return(cs.raw)
}


############################################################
# get the voxels that have any people with hemorrhage there
############################################################

ms = mat[which(rs > 0), ]
dist.mat = t(ms)
dist.mat = dist.mat %*% ms
A = matrix(diag(dist.mat), ncol=ncol(dist.mat), nrow=nrow(dist.mat))
At = t(A)

dice = 2* dist.mat / (At + A)
stopifnot(all(diag(dice) == 1))



device = "png"


lr_symm = function(img){
	dimg = dim(img)
	max.slice = dimg[1]	
	mid.slice = (max.slice+1)/2

	w = which(img > 0, arr.ind=TRUE)
  	## 20 - then 160, 90 - 20 + 90
  	## if 160 then 90 - 160 + 90
	w[, 1] = 2 * mid.slice - w[,1]
	w = w[ w[, 1] > 0 & w[, 1] < max.slice, ]
	img[w] = 1

	img = (img > 0)*1
	img = newnii(img)
}


nkeeps = c(30, 100, 500, 1000, 2000, 3000, c(0.05, .01, .001))



whichdir = "reoriented"
outcome = "NIHSS"
rerun = TRUE

for (outcome in c("NIHSS", "GCS")){
	adder = paste0(outcome, "_")
	if (outcome == "NIHSS"){
		adder = ""
	}




	outfile = file.path(outdir, paste0(adder, "Pvalue_Matrix.Rda"))
	load(file=outfile)

	### get top N voxels
	y = res[,"X","mod.1"]
	rrn = which(rs > ncut)
	rrn = rrn[order(y)]
	yord = y[order(y)]


	print(outcome)

	if (outcome == "NIHSS"){
		nkeep = .01
		addtext = "A"

	}
	if (outcome == "GCS"){
		nkeep = 1000
		addtext = "B"
	}
		orig = nkeep
	# nkeep = 3000
		if (orig < 1){
			nkeep = sum(yord <= orig)
		}

		rn = rrn[seq(nkeep)]


		outstub = file.path(outdir, 
			paste0(adder, "Top_", orig, "_pvalues_Final"))
		fp = paste0(outstub, ".", device)

		outimg = paste0(outstub, ".nii.gz")	
		fp_symm = paste0(outstub, "_symm.", device)
		outimg_symm = paste0(outstub, "_symm.nii.gz")	

		if (!file.exists(fp) | !file.exists(outimg) | rerun |
			!file.exists(outimg_symm) | !file.exists(fp_symm) ){
			open_dev(fp)
			res.p = temp
			res.p[!is.na(res.p)] = NA
			res.p[rn] = 1


			xyz= NULL
			xyz = floor(colMeans(which(res.p > 0, arr.ind = TRUE)))

			# res.p[ rs > ncut ]  = 1-y
			# cols = c("blue", "green", "yellow", "orange", "red", "white")
			mask.overlay(temp, res.p, col.y="red", 
				xyz= xyz, text=addtext, text.cex = 8,
					text.x= 32, text.y=32)
			dev.off()

			res.p[is.na(res.p)] = 0
			res.p = (res.p > 0)*1
			res.p = newnii(res.p)
			writeNIfTI(res.p, file=outstub)

			res.psymm = lr_symm(res.p)
			open_dev(fp_symm)
			print(addtext)
			mask.overlay(temp, res.psymm, col.y="red", 
				xyz= xyz, text=addtext, text.cex = 8,
					text.x= 32, text.y=32)
			dev.off()		
			writeNIfTI(res.psymm, file=paste0(outstub, "_symm"))
		}

}
