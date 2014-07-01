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
  if (!keepall) cs.raw = cs.raw[cs.raw != 0, , drop=FALSE]
  rownames(cs.raw) = names(ind.list)
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

nkeeps = c(30, 100, 500, 1000, 2000, 3000)

for (pval in c(0.05, .01, .001)){

	nkeeps = c(nkeeps, sum(yord <= pval))
}


for (nkeep in nkeeps){
	rn = rrn[seq(nkeep)]

	outstub = file.path(outdir, 
		paste0(adder, "Top_", nkeep, "_pvalues"))
	fp = paste0(outstub, ".", device)

	outimg = paste0(outstub, ".nii.gz")	

	if (!file.exists(fp) | !file.exists(outimg) | rerun){
		open.dev(fp)
		res.p = temp
		res.p[!is.na(res.p)] = NA
		res.p[rn] = 1


		xyz= NULL
		xyz = floor(colMeans(which(res.p > 0, arr.ind = TRUE)))

		# res.p[ rs > ncut ]  = 1-y
		# cols = c("blue", "green", "yellow", "orange", "red", "white")
		mask.overlay(temp, res.p, col.y="red", xyz= xyz, text=paste0("ROI for ", outcome, " score\n", "(V=", nkeep, ")"))
		dev.off()

		res.p[is.na(res.p)] = 0
		res.p = (res.p > 0)*1
		res.p = newnii(res.p)
		writeNIfTI(res.p, file=outstub)
	}
}




atfile = file.path(atlasdir, "All_FSL_Atlas_Labels.Rda")

load(file=atfile)
atfile = file.path(atlasdir, 
  paste0(whichdir, "_All_Atlas_ROI_Overlap_Measures.Rda"))
load(file=atfile)


lists = list(mni.list, jhut1.list, jhut2.list)

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

pop.pcts = sapply(sublists, function(indlist){
	sapply(indlist, function(x) mean(mat[x, ]))
})
colnames(pop.pcts) = c("EVE_1", "EVE_2")


# allres = allres
make.pvalimg = function(pvalimg){
	pvalimg.tab = llply(lists, function(x) {
		x = area_pct(pvalimg, ind.list=x, keepall=TRUE)		
		x$nvox = x$nvox/sum(x$nvox) * 100
		x$roi_pct_any = x$roi_pct_any * 100
		x$roi_mean_pct = x$roi_mean_pct * 100		
		x = x[order(x$nvox, decreasing=TRUE), , drop=FALSE]
		x$area = rownames(x)
		x
	}, .progress= "text")

	names(pvalimg.tab) = c("MNI", "EVE_1", "EVE_2")
	return(pvalimg.tab)
}


for (nkeep in c(1000, 2000, 3000)){

# nkeep = 3000
	rn = rrn[seq(nkeep)]
	pval = yord[nkeep]

	ppcts = sapply(sublists, function(indlist){
		sapply(indlist, function(x) {
			i = intersect(x, rn)
			length(i)/length(x)
		})
	})
	colnames(ppcts) = c("EVE_1", "EVE_2")


	submat = mat[rn,]
	wi = colSums(submat)

	pvalimg = array(0, dim = dim(jhut1.img))
	pvalimg[rn] = 1

	pvalimg.tab = make.pvalimg(pvalimg)



	outfile = file.path(outdir, 
		paste0(adder, "Top_", nkeep, "_Pvalues_df.Rda"))
	save(submat, rs, rn, wi, nkeep, 
		dist.mat, dice, pval, pvalimg.tab, ppcts, pop.pcts,
		file=outfile)
}

for (pval in c(0.05, .01, .001)){

	# nkeep = 3000
	nkeep = sum(yord <= pval)
	rn = rrn[seq(nkeep)]

	ppcts = sapply(sublists, function(indlist){
		sapply(indlist, function(x) {
			i = intersect(x, rn)
			length(i)/length(x)
		})
	})
	colnames(ppcts) = c("EVE_1", "EVE_2")


	submat = mat[rn,]
	wi = colSums(submat)

	pvalimg = array(0, dim = dim(jhut1.img))
	pvalimg[rn] = 1

	pvalimg.tab = make.pvalimg(pvalimg)


	outfile = file.path(outdir, 
		paste0(adder, "Top_", pval, "_Pvalues_df.Rda"))
	save(submat, rs, rn, wi, pval, nkeep,
		dist.mat, dice, pvalimg.tab, ppcts, pop.pcts,
		file=outfile)
}




# Top 2000 %
# Top 3000 %
# 0.05 0.01