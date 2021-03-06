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


outcome = "GCS"
adder = paste0(outcome, "_")
if (outcome == "NIHSS"){
	adder = ""
}

x = load(file.path(outdir, 
	paste0("Regress_ROI_", outcome, "_Results.Rda")))


template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)

t.t1 = file.path(tempdir, "sct1_unsmooth.nii.gz")
temp.t1 = readNIfTI(t.t1)

setwd(progdir)

ids.111 = read.csv(file.path(basedir, "111_patients.csv"), 
	stringsAsFactors= FALSE)
uid = ids.111$patientName
all.ids = ids.111$id


outfile = file.path(atlasdir, 
  paste0(whichdir, "_All_Atlas_ROI_Overlap_Measures.Rda"))

load(file=outfile)


nkeeps = vox.nkeeps
nkeep = nkeeps[1]

all.res = vector("list", length=length(nkeeps))
names(all.res) = paste0("vox.", nkeeps)

for (nkeep in nkeeps){
	outimg = file.path(outdir, 
		paste0(adder, "Top_", nkeep, "_pvalues.nii.gz"))

	popimg = readNIfTI(outimg)

	makeres = function(allres){
		allres$fname = gsub("^bws", "", allres$fname)
		allres$fname = gsub("ROI$", "", allres$fname)
		allres[, c("raw", "bin", "weighted")] = 
			round(allres[, c("raw", "bin", "weighted")] *100, 1)
		ss = strsplit(allres$fname, "_")
		allres$id = sapply(ss, function(x) x[1])
		return(allres)
	}

	names(mni.list)[names(mni.list) == "Uncategorized"] = "Ventricles/WM"
	names(jhut1.list)[names(jhut1.list) == "Background"] = "Ventricles"
	names(jhut2.list)[names(jhut2.list) == "Background"] = "Ventricles"


	# view.png(spm1_t1_hot_overlay.png)
	lists = list(mni.list, jhut1.list, jhut2.list)
	pop.tab = llply(lists, function(x) {
		tt = area_pct(popimg, ind.list=x, keepall=TRUE)
	}, .progress= "text")
	sums = sapply(pop.tab, function(x) sum(x$nvox))
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
	lname = paste0("vox.", nkeep)
	all.res[[lname]] = list(dfs=dfs, tops=tops, xtabs = xtabs)
}

save(all.res, nkeeps, outcome, file=file.path(outdir, 
	paste0(outcome, "_Engagement_Breakdown.Rda")))