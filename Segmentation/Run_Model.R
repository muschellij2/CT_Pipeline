#####################################################################
## This code is for Image Segmentation of CT
## The code is R but calls out FSL
##
## Author: John Muschelli
## Last updated: May 20, 2014
#####################################################################
#####################################################################
rm(list=ls())
library(plyr)
library(cttools)
library(ROCR)
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

segdir = file.path(progdir, "Segmentation")
source(file.path(segdir, "performance_functions.R"))

#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )

fnames = gsub("^bws", "", fnames)
fnames = paste0(fnames, ".gz")
ids = gsub("(\\d\\d\\d-(\\d|)\\d\\d\\d)_.*", "\\1", fnames)
fdf = data.frame(id = ids, stringsAsFactors= FALSE)
fdf$iddir = file.path(basedir, fdf$id)
fdf$outdir = file.path(fdf$iddir, "Predictors")
makedir = sapply( fdf$outdir, dir.create, showWarnings =FALSE)
fdf$roi = file.path(rootdir, "ROI_data", fdf$id, fnames)
fdf$img = file.path(fdf$iddir, gsub("ROI\\.nii", ".nii", fnames))
fdf$mask = file.path(fdf$iddir, 
	"Skull_Stripped", 
	gsub("ROI\\.nii", "_SS_Mask_0.01.nii", fnames))
irow = 2
x = fdf[irow,]



# fdf = fdf[1:10,]
fpr.stop = .1
run_model = function(x, fpr.stop = .1){

	system.time({
		img.pred = make_predictors(
		img=x$img, 
		mask = x$mask, 
		roi = x$roi,
		nvoxels = 1, 
		moments = 1:4, lthresh = 40, uthresh = 80,
		save_imgs = TRUE, 
		outdir = x$outdir,
		overwrite = FALSE, 
		verbose= TRUE)
	})
	df = img.pred$df
	keep.ind = img.pred$keep.ind
	nim = img.pred$nim
	df = df[ keep.ind, ]
	# df$Y = roi[keep.ind]

	# img = readNIfTI(x$img, reorient= FALSE)
	# # mask = readNIfTI(x$mask, reorient= FALSE)
	# # erode the mask
	# mask = fslerode(file=x$mask, kopts = "-kernel box 1x1x1", 
	#                   reorient=FALSE, retimg = TRUE)
	# roi = readNIfTI(x$roi, reorient = FALSE)
	# mask = mask > 0
	# mask[ roi == 1] = 1


	samps = seq(nrow(df)) %in% sample(nrow(df), size= 1e4)
	train = df[samps,]
	test = df[!samps,]

	mod = glm(Y ~ ., data=train, family=binomial())

	test.pred = predict(mod, newdata=test, type="response")
	# train.pred = predict(mod, newdata=train, type="response")

	pred <- prediction( test.pred, test$Y)
	# perf <- performance(pred,"tpr","fpr")
	# xind = perf@x.values[[1]] <= fpr.stop
	# perf@x.values[[1]] = perf@x.values[[1]][xind]
	# perf@y.values[[1]] = perf@y.values[[1]][xind]

	# auc = performance(pred, "auc")@y.values[[1]]
	# plot(perf)
	pauc = performance(pred, "auc", fpr.stop= fpr.stop)
	pauc = pauc@y.values[[1]] / fpr.stop
	# pauc

	smod = summary(mod)
	smod$deviance.resid = NULL
	mod = remove_lmparts(mod)

	# rownames(df) = NULL

	return(list(mod=mod, pauc = pauc, 
		train.ind = samps, 
		smod = smod,
		keep.ind = keep.ind, nim= nim))
}

mods = alply(fdf, 1, 
	run_model, 
	.progress = "text")

save(mods, file = file.path(outdir, "Segmentation_Models.Rda"))


# plot(perf)







# # function(x){
# irow = 5
# x = fdf[irow,]
# new.img = readNIfTI(x$img, reorient= FALSE)
# new.mask = readNIfTI(x$mask, reorient= FALSE)
# # erode the mask
# new.mask = fslerode(file=new.mask, kopts = "-kernel box 1x1x1", 
#                   reorient=FALSE, retimg = TRUE)

# new.roi = readNIfTI(x$roi, reorient = FALSE)

# ind = which(new.roi > 0 & new.mask < 1, arr.ind = TRUE)

# new.img.pred = make_predictors(img=new.img, mask = new.mask, nvoxels = 1, moments = 1:4,
#                             lthresh = 40, uthresh = 80)
# new.df = new.img.pred
# new.df = data.frame(new.df)
# new.df$Y = c(new.roi)

# new.pred = predict(mod, newdata=new.df, type="response")
# pimg = new.roi
# pimg@.Data = array(new.pred, dim=dim(pimg))
# pimg = cal_img(pimg)
# pimg[new.mask < 1 ] = 0

# run.pred = new.pred[new.mask > 0 | new.df$Y > 0]
# Y = new.df$Y[new.mask > 0 | new.df$Y > 0]
# # centroid = which(mask > 0, arr.ind=TRUE)
# # dcent = t(t(centroid) - colMeans(centroid))
# # dcent = sqrt(rowSums(dcent^2))

# pred <- prediction( run.pred, Y)
# perf <- performance(pred,"tpr","fpr")


# auc = performance(pred, "auc")@y.values[[1]]
# # plot(perf)
# fpr.stop = .1
# pauc = performance(pred, "auc", fpr.stop= fpr.stop)
# pauc = pauc@y.values[[1]] / fpr.stop
# pauc



