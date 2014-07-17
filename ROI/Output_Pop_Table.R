rm(list=ls())
library(fslr)
library(xtable)
library(cttools)
### need cairo for cluster
options(bitmapType='cairo')
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

outdir = file.path(basedir, "results")
whichdir = "reoriented"

tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

load(file.path(atlasdir, "All_FSL_Atlas_Labels.Rda"))

# whichdir = "FLIRT"
rerun = FALSE

template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)

t.t1 = file.path(tempdir, "sct1_unsmooth.nii.gz")
temp.t1 = readNIfTI(t.t1)

setwd(progdir)

ids.111 = read.csv(file.path(basedir, "111_patients.csv"), 
	stringsAsFactors= FALSE)
uid = ids.111$patientName
all.ids = ids.111$id


spm.popimg = file.path(outdir, paste0(whichdir, "_Weighted_Sum_Image"))
spm.binimg = file.path(outdir, paste0(whichdir, "_Binary_Sum_Image"))

spm.popimg = paste0(spm.popimg, ".nii.gz")
spm.binimg = paste0(spm.binimg, ".nii.gz")

popimg = readNIfTI(spm.popimg)

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

names(mni.list)[names(mni.list) == "Uncategorized"] = "CSF/WM"
names(jhut1.list)[names(jhut1.list) == "Background"] = "CSF"
names(jhut2.list)[names(jhut2.list) == "Background"] = "CSF"


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



# view.png(spm1_t1_hot_overlay.png)
lists = list(mni.list, jhut1.list, jhut2.list)
names(lists) = c("MNI", "EVE_1", "EVE_2")
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


pop.tab = llply(lists, function(x) {
	tt = area_pct(popimg, ind.list=x, keepall=TRUE)
}, .progress= "text")
sums = sapply(pop.tab, function(x) sum(x$nvox))
stopifnot(all(diff(sums) == 0))

col.pop.tab = llply(col.lists, function(x) {
	tt = area_pct(popimg, ind.list=x, keepall=TRUE)
}, .progress= "text")
sums = sapply(col.pop.tab, function(x) sum(x$nvox))
stopifnot(all(diff(sums) == 0))

##### scaling them to %
pop.tab = llply(pop.tab, function(x) {
	x$nvox = x$nvox/sum(x$nvox) * 100
	x$roi_pct_any = x$roi_pct_any * 100
	x$roi_mean_pct = x$roi_mean_pct * 100
	x = x[order(x$nvox, decreasing=TRUE), , drop=FALSE]
	x$area = rownames(x)
	x
}, .progress= "text")

nm = names(pop.tab) = c("MNI", "EVE_1", "EVE_2")

##### scaling them to %
col.pop.tab = llply(col.pop.tab, function(x) {
	x$nvox = x$nvox/sum(x$nvox) * 100
	x$roi_pct_any = x$roi_pct_any * 100
	x$roi_mean_pct = x$roi_mean_pct * 100
	x = x[order(x$nvox, decreasing=TRUE), , drop=FALSE]
	x$area = rownames(x)
	x
}, .progress= "text")

col.nm = names(col.pop.tab) = c("EVE_1", "EVE_2")

dfs = tops = xtabs = 
	vector(mode = "list", length=length(pop.tab))
for (itab in seq(pop.tab)){
	df = allres = pop.tab[[itab]]
	df$nvox = round(df$nvox, 2)
	df$roi_mean_pct = round(df$roi_mean_pct, 2)
	df$roi_pct_any = round(df$roi_pct_any, 2)
	df = df[df$nvox != 0, , drop=FALSE]
	top.area = rownames(df)[1]
	xdf = df
	df = df[, c("area", "nvox")]		

	colnames(df) = c("Area", nm[itab])
	xtabs[[itab]] = xtable(df)
	tops[[itab]] = top.area 
	xdf = xdf[, c("area", "nvox", "roi_pct_any", "roi_mean_pct")]		
	colnames(xdf) = c("Area", nm[itab], 
		paste0(nm[itab], "_ROI_Pct_Any"),
		paste0(nm[itab], "_ROI_Pct"))
	dfs[[itab]] = xdf
}
names(dfs) = names(tops) = names(xtabs) = nm
id = "Population"


col.dfs = col.tops = col.xtabs = 
	vector(mode = "list", length=length(col.pop.tab))
for (itab in seq(col.pop.tab)){
	df = allres = col.pop.tab[[itab]]
	df$nvox = round(df$nvox, 2)
	df$roi_mean_pct = round(df$roi_mean_pct, 2)
	df$roi_pct_any = round(df$roi_pct_any, 2)
	df = df[df$nvox != 0, , drop=FALSE]
	top.area = rownames(df)[1]
	xdf = df
	df = df[, c("area", "nvox")]		

	colnames(df) = c("Area", nm[itab])
	col.xtabs[[itab]] = xtable(df)
	col.tops[[itab]] = top.area 
	xdf = xdf[, c("area", "nvox", "roi_pct_any", "roi_mean_pct")]		
	colnames(xdf) = c("Area", nm[itab], 
		paste0(nm[itab], "_ROI_Pct_Any"),
		paste0(nm[itab], "_ROI_Pct"))
	col.dfs[[itab]] = xdf
}
names(col.dfs) = names(col.tops) = names(col.xtabs) = col.nm
id = "Population"


save(dfs, xtabs, tops, col.dfs, col.xtabs, col.tops, 
	file=file.path(outdir, "Population_Table.Rda"))