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

whichdir = "reoriented"
outcome = "NIHSS"
adder = paste0(outcome, "_")
if (outcome == "NIHSS"){
	adder = ""
}
stopifnot(outcome == "NIHSS")


rerun = FALSE


outfile = file.path(outdir, paste0(adder, "Pvalue_Matrix.Rda"))
load(file=outfile)

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

### get top N voxels
y = res[,"X","mod.1"]
rrn = which(rs > ncut)
rrn = rrn[order(y)]
yord = y[order(y)]





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


# for (pval in c(0.05, .01, .001)){

# 	nkeeps = c(nkeeps, sum(yord <= pval))
# }

nkeeps = 19047
nkeep = nkeeps[1]
for (nkeep in nkeeps){

	orig = nkeep
# nkeep = 3000
	if (orig < 1){
		nkeep = sum(yord <= orig)
	}

	rn = rrn[seq(nkeep)]


	rn = rrn[seq(nkeep)]

	outstub = file.path(outdir, 
		paste0(adder, "Top_", orig, "_pvalues"))
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
			xyz= xyz, text=paste0("ROI for ", outcome, 
				" score\n", "(V=", nkeep, ")"))
		dev.off()

		res.p[is.na(res.p)] = 0
		res.p = (res.p > 0)*1
		res.p = newnii(res.p)
		writeNIfTI(res.p, file=outstub)

		res.psymm = lr_symm(res.p)
		open_dev(fp_symm)

		mask.overlay(temp, res.psymm, col.y="red", 
			xyz= xyz, 
			text=paste0("Symmetrized ROI for ", outcome, 
				" score\n", "(V=", nkeep, ")"))
		dev.off()		
		writeNIfTI(res.psymm, file=paste0(outstub, "_symm"))
	}
}




atfile = file.path(atlasdir, "All_FSL_Atlas_Labels.Rda")

load(file=atfile)
atfile = file.path(atlasdir, 
  paste0(whichdir, "_All_Atlas_ROI_Overlap_Measures.Rda"))
load(file=atfile)


lists = list(mni.list, jhut1.list, jhut2.list)
names(lists) = c("MNI", "EVE_1", "EVE_2")

sublists = list(jhut1.list, jhut2.list)

sublists = lapply(sublists, function(x) {
	area = names(x)
	x[grep("GLOBUS_PALLIDUS|THALAMUS|PUTAMEN", area)]
})

sublists = lapply(sublists, function(x) {
	xx = unlist(x)
	area = names(xx)
	area = gsub("_left\\d*", "", area)
	area = gsub("_right\\d*", "", area)
	uarea = unique(area)
	x = lapply(uarea, function(aname){
		ind = which(area %in% aname)
		xx[ind]
	})
	names(x) = uarea
	x
})

names(sublists) = c("EVE_1", "EVE_2")

pop.pcts = sapply(sublists, function(indlist){
	sapply(indlist, function(x) mean(mat[x, ]))
})
colnames(pop.pcts) = c("EVE_1", "EVE_2")


col.lists = list(jhut1.list, jhut2.list)
names(col.lists) = c("EVE_1", "EVE_2")

col.lists = lapply(col.lists, function(x) {
	area = names(x)
	area = gsub("_left", "", area)
	area = gsub("_right", "", area)
	uarea = unique(area)
	res = lapply(uarea, function(aname){
		ind = which(area %in% aname)
		xx = sort(unlist(x[ind]))
		names(xx) = NULL
		print(ind)
		xx
		# xx[ind]
	})
	names(res) = uarea
	res
})

pop.colpcts = sapply(col.lists, function(indlist){
	sapply(indlist, function(x) mean(mat[x, ]))
})
names(pop.colpcts) = c("EVE_1", "EVE_2")


# allres = allres
make.pvalimg = function(pvalimg, runlist = lists){
	pvalimg.tab = llply(runlist, function(x) {
		x = area_pct(pvalimg, ind.list=x, keepall=TRUE)		
		x$nvox = x$nvox/sum(x$nvox) * 100
		x$roi_pct_any = x$roi_pct_any * 100
		x$roi_mean_pct = x$roi_mean_pct * 100		
		x = x[order(x$nvox, decreasing=TRUE), , drop=FALSE]
		x$area = rownames(x)
		x
	}, .progress= "text")

	names(pvalimg.tab) = names(runlist)
	return(pvalimg.tab)
}


for (nkeep in nkeeps){

	orig = nkeep
# nkeep = 3000
	if (orig < 1){
		nkeep = sum(yord <= orig)
	}
	rn = rrn[seq(nkeep)]

	#### get indices for symmetrized version
	res.p = temp
	res.p[!is.na(res.p)] = NA
	res.p[rn] = 1

	res.p[is.na(res.p)] = 0
	res.p = (res.p > 0)*1

	res.psymm = lr_symm(res.p)
	rn_symm = which(res.psymm > 0)

	if (orig < 1){
		pval = orig
	} else {
		pval = yord[nkeep]
	}

	ppcts = sapply(sublists, function(indlist){
		sapply(indlist, function(x) {
			i = intersect(x, rn)
			length(i)/length(x)
		})
	})
	colnames(ppcts) = c("EVE_1", "EVE_2")

	col.ppcts = sapply(col.lists, function(indlist){
		sapply(indlist, function(x) {
			i = intersect(x, rn)
			length(i)/length(x)
		})		
	})
	names(col.ppcts) = c("EVE_1", "EVE_2")


	ppcts_symm = sapply(sublists, function(indlist){
		sapply(indlist, function(x) {
			i = intersect(x, rn_symm)
			length(i)/length(x)
		})
	})
	colnames(ppcts_symm) = c("EVE_1", "EVE_2")


	submat = mat[rn,]
	wi = colSums(submat)

	pvalimg = array(0, dim = dim(jhut1.img))
	pvalimg[rn] = 1

	pvalimg.tab = make.pvalimg(pvalimg, lists)

	col.pvalimg.tab = make.pvalimg(pvalimg, col.lists)

	pvalimg = array(0, dim = dim(jhut1.img))
	pvalimg[rn_symm] = 1

	pvalimg.tab_symm = make.pvalimg(pvalimg, lists)

	submat_symm = mat[rn_symm,]
	wi_symm = colSums(submat_symm)

	outfile = file.path(outdir, 
		paste0(adder, "Top_", orig, "_Pvalues_df.Rda"))
	save(submat, rs, rn, wi, nkeep, 
		dist.mat, dice, pval, 
		pvalimg.tab, pvalimg.tab_symm, col.pvalimg.tab, 
		ppcts, pop.pcts, col.ppcts,
		ppcts_symm,
		wi_symm, submat_symm,
		file=outfile)
}
