#################################
# Symmetric Analysis
#################################
rm(list=ls())
library(cttools)
library(oro.nifti)
library(scales)
library(RColorBrewer)
library(data.table)
library(ggplot2)
library(grid)
library(plyr)
library(smallPDF)
homedir = "/Applications"
rootdir = "~/CT_Registration"
basedir = file.path(rootdir, "data")
outdir = basedir
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
  basedir = file.path(rootdir, "Registration")
  outdir = file.path(basedir, "results")
}
progdir = file.path(rootdir, "programs")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")


get.id = function(x){
  ss = strsplit(x, "_")
  ss = sapply(ss, head, 1)
  ss = gsub(".*(\\d\\d\\d-.*)", "\\1", ss)
  ss
}

id_to_pname = function(x){
  as.numeric(gsub("-", "", x))
}

demog = read.csv(file=file.path(basedir, "Demog_NIHSS_Mask.csv"), 
 stringsAsFactors=FALSE)
demog$Base_ICH_10 = demog$Diagnostic_ICH /10


template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)
dtemp = dim(temp)

t1 = file.path(tempdir, "betsct1_unsmooth.nii.gz")
temp.t1 = readNIfTI(t1)

#### load voxel data
outfile = file.path(outdir, "Voxel_Matrix.Rda")
vmat = load(file=outfile )

stopifnot(nrow(mat) == prod(dtemp))

nim = temp
nim[is.na(nim)] = 0
nim[!is.na(nim)] = 1

max.slice = dtemp[1] 
mid.slice = (max.slice+1)/2

w = which(nim > 0, arr.ind=TRUE)
ind = which(nim > 0)
w = cbind(w)
w = data.frame(w)


w2 = w
	## 20 - then 160, 90 - 20 + 90
	## if 160 then 90 - 160 + 90
w2[, 1] = 2 * mid.slice - w2[,1]
w2 = w2[ w2[, 1] > 0 & w2[, 1] <= max.slice, ]

colnames(w2) = paste0("flip.", colnames(w2))
w$ind = ind



stopifnot(nrow(mat) == nrow(w))

##### let's choose the left side 

w$left = w$dim1 > mid.slice
w$left[ w$dim1 == mid.slice ] = NA

ww = cbind(w, w2)

mids = ww[ ww$dim1 == ww$flip.dim1 & ww$dim2 == ww$flip.dim2 &
	ww$dim3 == ww$flip.dim3, ]

left = ww[ which(ww$left == TRUE), ]

dat = rbind(mids, left)

dat$left = NULL

colnames(dat) = gsub("^dim", "orig.dim", colnames(dat))
colnames(dat) = gsub("^flip\\.", "", colnames(dat))

dat = merge(dat, w, all.x=TRUE, sort=FALSE, by=paste0("dim", 1:3), 
	suffixes = c(".flip",".orig"))

dat = dat[ order(dat$dim1, dat$dim2, dat$dim3), ]
rownames(dat) = NULL

dat$left = NULL


symm_dat = mat[ dat$ind.orig, ]
flip_dat = mat[ dat$ind.flip, ]

stopifnot(all.equal(dim(symm_dat), dim(flip_dat)))

symm = symm_dat | flip_dat

dat$symm_rs = rowSums(symm)
cs = colSums(mat)
symm_cs = colSums(symm)
ncut = 10

keep = dat$symm_rs > ncut
symm_mat = symm[ keep, ]

nihss.mods = fast_lm(Y = demog$Enrollment_NIHSS_Total, 
	X = t(symm_mat), ncheck = 1e3, spot.check=TRUE)

gcs.mods = fast_lm(Y = demog$Enrollment_GCS_Add, 
	X = t(symm_mat), ncheck = 1e3, spot.check=TRUE)


mods = list(NIHSS=nihss.mods, GCS=gcs.mods)

results = dat[ keep, ]

# for (pval in c(0.05, .01, .001)){

# 	nkeeps = c(nkeeps, sum(yord <= pval))
# }

results = cbind(results, NIHSS=nihss.mods$p.val,
	GCS = gcs.mods$p.val)

all.rn = dat$ind.orig[dat$symm_rs > 0]

nkeeps = c(1000, 2000, 3000, c(0.05, .01, .001), 
	nrow(results), length(all.rn))
nkeep = nrow(results)

outcome = "NIHSS"

#### To get symm picture
res.p = res.orig = temp

# for (i in 80:85){
# 	res.p[!is.na(res.p)] = NA
# 	res.p[ dat$ind.orig] = symm[,i]
# 	res.p[ res.p == 0] = NA
# 	# ortho2(temp, res.p)

# 	res.orig[!is.na(res.orig)] = NA
# 	res.orig[seq(prod(dim(res.orig)))] =  mat[,i]
# 	res.orig[ res.orig == 0] = NA
# 	ortho2(res.p, res.orig, col.y="red")
# }


for (outcome in c("NIHSS", 'GCS')){
	y = results[, outcome]
	ord = order(y)
	res = results[ord, ]
	rrn = res$ind.orig
	yord = y[ord]

	for (nkeep in nkeeps){

		orig = nkeep
	# nkeep = 3000
		if (orig < 1){
			nkeep = sum(yord <= orig)
		}
		if (nkeep == length(all.rn)){
			rn = all.rn
		} else { 
			rn = rrn[seq(nkeep)]
		}


		outstub = file.path(outdir, 
			paste0(outcome, "_Symm_analysis_Top_", orig, "_pvalues"))
		fp = paste0(outstub, ".png")

		# outimg = paste0(outstub, ".nii.gz")	
		# fp_symm = paste0(outstub, "_symm.", device)
		# outimg_symm = paste0(outstub, "_symm.nii.gz")	

		# if (!file.exists(fp) | !file.exists(outimg) | rerun |
		# 	!file.exists(outimg_symm) | !file.exists(fp_symm) ){
			open_dev(fp)
			res.p = temp
			res.p[!is.na(res.p)] = NA
			res.p[rn] = 1
			# res.p[rn] = yord[seq(nkeep)]


			xyz= NULL
			xyz = floor(colMeans(which(res.p > 0, arr.ind = TRUE)))

			# res.p[ rs > ncut ]  = 1-y
			# cols = c("blue", "green", "yellow", "orange", "red", "white")
			mask.overlay(temp, res.p, col.y="red", 
				xyz= xyz, text=paste0("ROI for ", outcome, 
					" score\n", "(V=", nkeep, ")"))
			dev.off()
			if ( nkeep == nrow(results)) {
				outstub = file.path(outdir, 
					paste0(outcome, "_Symm_analysis_Top_", orig, 
						"_pvalue_map"))
				fp = paste0(outstub, ".png")


				thresholds = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1)
				res.p = temp
				res.p[!is.na(res.p)] = NA
				stopifnot(all(y > 0))
				cuts = cut(yord, breaks=thresholds, 
					include.lowest=FALSE)
				cuts = factor(cuts, levels=rev(levels(cuts)))
				levs = levels(cuts)
				cuts = as.numeric(cuts)
				stopifnot(length(rrn) == length(cuts))
				res.p[ rrn ] = cuts
				# res.p[rs > ncut] = -log(y)
				# res.p[rs > ncut] = as.numeric(y <= 0.05)
				# res.p[ rs > ncut ]  = 1-y
				cols = c("blue", "deepskyblue1", "green", "yellow", 
					"orange", "red")

				##### all this depends on what unique cuts are shown in the
				#### perspective
				xyz = c(91, 109, 91)
	

				cols = alpha(cols, .7)

				open_dev(fp)

				mask.overlay(temp.t1, res.p, window = c(0, 1000),
					col.y=cols, zlim.y = c(0, 6),
					addlegend = TRUE,
					leg.x = 15, leg.y= 60, 
					legend=rev(levs), 
					leg.col=alpha(rev(cols), 1), leg.cex=1.5,
					leg.title = "P-value"
				)

				# res.p[rn] = yord[seq(nkeep)]

				# ortho2(temp, res.p, xyz=xyz)
				dev.off()

				adj.p = p.adjust(yord, method="bonferroni")
				sig.inds = adj.p < .05
				if (any(sig.inds)) {
					print(paste0(outcome, ": Significance!"))
					print(paste0(sum(sig.inds), " Significant Voxels"))
				} else {
					print(paste0(outcome, ": No Significance"))
				}
				# y[]
			}
		print(nkeep)
	}
}

# Bonf Sig for NIHSS
# res.p = temp
# res.p[!is.na(res.p)] = NA
# ind = c(2731540L, 3714373L, 2848296L, 2810638L, 2849915L, 3674916L, 
# 2927744L, 2967564L, 3006841L, 2849553L, 2888830L, 2928107L, 2967384L, 
# 3006661L, 2769741L, 2809018L, 3477624L, 3163588L, 2810637L, 3006480L, 
# 2809019L, 3753469L, 2849733L, 2889010L, 3045757L, 3163769L, 3045758L, 
# 3085035L, 2887573L, 2888829L, 3006481L, 2849735L, 3045938L, 2888467L, 
# 2928106L, 2967383L, 3006660L, 3714192L, 2967567L, 2928290L, 2928287L, 
# 2849914L, 2889192L, 2928469L, 2967746L, 3007023L, 2809913L, 2771541L, 
# 2850096L, 2888648L, 2927925L, 2771179L, 2810456L, 3045939L, 3085216L, 
# 3714554L, 3675097L, 2770817L, 2810094L, 3046118L, 3085395L, 2692263L, 
# 3085034L, 3124311L, 3124492L, 2810457L, 2849734L, 2889012L, 3006662L, 
# 2889191L, 2928468L, 2967745L, 3007022L, 2887391L, 2926668L, 2965945L, 
# 2887572L, 2926849L, 3557808L, 2810818L, 2889011L, 2928288L, 2967565L, 
# 2965944L, 2888651L, 2888650L, 2927927L, 2967204L, 2928108L, 2967385L, 
# 2850638L, 2889735L, 2769561L, 2849190L, 3085215L, 2967386L, 3006663L, 
# 3006844L, 3516901L, 2847752L, 3714553L, 2770998L, 2810275L, 2889373L, 
# 3753830L, 2929193L, 2968470L, 2849374L, 2850277L, 2889554L, 2811361L, 
# 2888106L, 2810095L, 2810276L, 2967205L, 3438528L, 2849554L, 2770457L
# )
# res.p[ind] = 1
# res.p[is.na(res.p)] = 0

# whichdir = "reoriented"
# outfile = file.path(atlasdir, 
#   paste0(whichdir, "_All_Atlas_ROI_Overlap_Measures.Rda"))

# load(file=outfile)

# load(file.path(atlasdir, "All_FSL_Atlas_Labels.Rda"))

# names(mni.list)[names(mni.list) == "Uncategorized"] = "CSF/WM"
# names(jhut1.list)[names(jhut1.list) == "Background"] = "CSF"
# names(jhut2.list)[names(jhut2.list) == "Background"] = "CSF"


# area_pct = function(img, ind.list, keepall) {
#   ## get overlap of indices
#   raw.mat = sapply(ind.list, function(x) sum(img[x]))
#   any.mat = sapply(ind.list, function(x) mean(img[x] > 0))
#   mn.mat = sapply(ind.list, function(x) mean(img[x]))
#   names(raw.mat) = names(ind.list)
#   ## cs is sum of indices of overlap
#   cs.raw = data.frame(nvox=raw.mat, roi_pct_any = any.mat,
#   	roi_mean_pct = mn.mat) 
#   rownames(cs.raw) = names(ind.list)
#   if (!keepall) cs.raw = cs.raw[rowSums(cs.raw) != 0, , drop=FALSE]
#   return(cs.raw)
# }



# # view.png(spm1_t1_hot_overlay.png)
# lists = list(mni.list, jhut1.list, jhut2.list)
# names(lists) = c("MNI", "EVE_1", "EVE_2")

# res.tab = llply(lists, function(x) {
# 	tt = area_pct(res.p, ind.list=x, keepall=FALSE)
# }, .progress= "text")
# sums = sapply(res.tab, function(x) sum(x$nvox))
# stopifnot(all(diff(sums) == 0))


# ##### scaling them to %
# res.tab = llply(res.tab, function(x) {
# 	x$nvox = x$nvox/sum(x$nvox) * 100
# 	x$roi_pct_any = x$roi_pct_any * 100
# 	x$roi_mean_pct = x$roi_mean_pct * 100
# 	x = x[order(x$nvox, decreasing=TRUE), , drop=FALSE]
# 	x$area = rownames(x)
# 	x
# }, .progress= "text")

# nm = names(res.tab) = c("MNI", "EVE_1", "EVE_2")