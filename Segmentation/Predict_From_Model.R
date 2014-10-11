###################################################################
## This code is for prediciton of Image Segmentation of CT
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

correct = "N3"
options = c("none", "N3", "N4", "N3_SS", "N4_SS",
        "SyN", "SyN_sinc", "Rigid", "Affine")
keep.obj = ls()

for (correct in options){
	all.obj = ls()
	rm.obj = all.obj[!(all.obj %in% c(keep.obj, "keep.obj"))]
	rm(list=rm.obj)
    for (i in 1:3) gc()
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

	short_predict = function(object, newdata, 
		lthresh=  .Machine$double.eps^0.5){
		tt <- terms(object)
		Terms <- delete.response(tt)
	    m <- model.frame(Terms, newdata, 
	    	na.action = na.pass, 
	        xlev = object$xlevels)
	    if (is.null(cl <- attr(Terms, "dataClasses"))) 
	            stop("no dataclasses")
	    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
	    # p <- object$rank
	    beta <- object$coefficients
	    beta = beta[ !is.na(beta) ]
	    predictor = drop(X[, names(beta), drop=FALSE ] %*% beta)
	   
		predictor <- family(object)$linkinv(predictor)
		predictor[ predictor < lthresh] = 0
		predictor
	}

	#### load voxel data
	outfile = file.path(outdir, "Voxel_Info.Rda")
	load(file=outfile )

	outfile = file.path(outdir, "111_Filenames.Rda")
	load(file = outfile)


	makedir = sapply( fdf$outdir, function(x) {
		if (!file.exists(x)){
			dir.create(x, showWarnings =FALSE)
		}
	})
	irow = 1
	x = fdf[irow,]

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
	mod.filename = file.path(outdir, 
		paste0("Collapsed_Models", adder, ".Rda"))
	load(mod.filename)
	vol.sdatas = vol.datas = vol.data
	reses = sreses = res

	get.pred <- as.numeric(Sys.getenv("SGE_TASK_ID"))
	if (is.na(get.pred)) get.pred = 112
	x = fdf[get.pred,]


# for (get.pred in runpreds){

	iddir = fdf$iddir[get.pred]
	id.outdir = fdf$outdir[get.pred]
	predname = nii.stub(basename(fdf$img[get.pred]))
	predname = file.path(id.outdir, 
		paste0(predname, "_predictors", adder, ".Rda"))
	load(predname)
	df = img.pred$df
	nim = img.pred$nim
	keep.ind = img.pred$keep.ind
	rm(list="img.pred")
    for (i in 1:3) gc()	
	df$in0100 = df$value >= 0 & df$value <= 100
	df$mask = df$mask > 0
	
	######################################
	# Get volume of ROI
	######################################	
	vdim = voxdim(nim)
	vres = prod(vdim) / 1000

	Y =  df$Y
	vol.roi = sum(Y) * vres	

	not0100 = sum(df$Y[ !df$in0100 ])

	######################################
	# Keep all ROI = 1, even if not inmask
	######################################	
	roi.not.in = which(df$Y == 1)
	roi.not.in = roi.not.in[!(roi.not.in %in% keep.ind)]
	keep.ind = sort(c(keep.ind, roi.not.in))


	#### need this because the length of df has changed
	roi.not.in = which(keep.ind %in% roi.not.in)	

	df = df[keep.ind,]
	Y =  df$Y
	benchmark = 1-mean(df$Y)
	brain.vol = nrow(df) * vres


	df$img = fdf$img[1]

	# for (imod in  seq(lmod)){
	# 	irow = which(vol.diff$imod == imod & 
	# 		vol.diff$get.pred == get.pred)

		preds = t(laply(all.mods, short_predict, newdata= df, 
                  .progress = "text"))

		#### we only think values 0 to 100 are actually blood
		rowMean = rowMeans(preds)
		rowMed = rowMedians(preds)
		rowMax = rowMaxs(preds)
		rowMin = rowMins(preds)
		rowProd = exp(rowSums(log(preds)))
		rowGeom = exp(rowMeans(log(preds)))

		colnames(preds) = colnames(res)[seq(ncol(preds))]
		preds = cbind(preds, rowMean, rowMed, rowMax, rowMin,
			rowProd, rowGeom)
		rownames(preds) = NULL

		#### if the mask did not contain these voxels, 
		#### then they are 0
		preds[roi.not.in, ] = 0

		#### only values 0 to 100 are likely to be ROI
		preds = preds * df$in0100


		vols = colSums(preds)

		vols = vols * vres
		vol.data[get.pred, ] = c(vol.roi, vols)

		# pred.imgs = vector(mode="list", length=ncol(preds))
		spreds = matrix(nrow=nrow(preds), ncol=ncol(preds))
		pb = txtProgressBar(max=ncol(preds), style=3)
		for (ipred in seq(ncol(preds))){
			nd = preds[, ipred]
			img = nim
			nd[roi.not.in] = 0
			img@.Data[keep.ind] = nd
			img = cal_img(img)
			img[is.na(img)]= 0
			img = datatype(img, 
				datatype= convert.datatype()$FLOAT32, 
				bitpix = convert.bitpix()$FLOAT32)
			cn = colnames(preds)[ipred]
			# pred.imgs[[ipred]] = img
			outimg = nii.stub(fdf$img[get.pred], bn=TRUE)
			outimg = file.path(id.outdir, 
				paste0(outimg, "_", cn, adder))			
			writeNIfTI(img, filename = outimg )

			#### smooth the results
			sm.img  = mean_image(img, nvoxels = 1)
			# sm.img = local_moment(img, nvoxels = 1, moment =1,
			# 	center = FALSE)[[1]]
			sm.img[abs(sm.img) <  .Machine$double.eps^0.5 ] = 0
			sm.img = sm.img[keep.ind]
			spreds[, ipred] = sm.img
			spreds[roi.not.in, ipred] = 0

			sm.img = nim
			sm.img@.Data[keep.ind] = spreds[, ipred]
			sm.img[is.na(sm.img)]= 0
			sm.img[roi.not.in] = 0

			sm.img = cal_img(sm.img)
			sm.img  = datatype(sm.img , 
				datatype= convert.datatype()$FLOAT32, 
				bitpix = convert.bitpix()$FLOAT32)

			outimg = paste0(outimg, "_smoothed", adder)
			writeNIfTI(sm.img, filename = outimg )
			setTxtProgressBar(pb, ipred)

		}
		close(pb)

		# spreds = laply(pred.imgs, function(img){
		# 	sm.img  = mean_image(img, 1)
		# 	sm.img[abs(sm.img) <  .Machine$double.eps^0.5 ] = 0
		# 	nd.smooth = sm.img[keep.ind]
		# 	return(nd.smooth)			
		# }, .progress = "text")
		# spreds = t(spreds)


		# rowMean = rowMeans(spreds)
		# rowMed = rowMedians(spreds)
		# spreds = cbind(spreds, rowMean, rowMed)

		##################
		# Deleted for space issues
		# predname = nii.stub(basename(fdf$img[get.pred]))
		# predname = file.path(outdir, 
		# 	paste0(predname, "_predictions", adder, ".Rda"))
		# save(preds, spreds, Y, file=predname, compress=TRUE)
		##################


		svols = colSums(spreds)
		svols = svols * vres

		vol.sdata[get.pred, ] = c(vol.roi, svols)


		# print(vol.roi)
		# print(vol.pred)		
		# print(vol.spred)		

		# roi = mods[[get.mod]]$nim
		# roi@.Data[keep.ind] = df$Y
		# roi = cal_img(roi)

		# mask.overlay(roi, img > .1, col.y="red", window=c(0, 1))

		fpr.stop = .1

		r.preds = alply(preds, 2, function(nd){
			pred <- prediction( nd, Y)			
		}, .progress = "text")

		accs = t(laply(r.preds, function(pred){
			acc = performance(pred, "acc")
			ind = which.max(acc@y.values[[1]])
			cutoff = acc@x.values[[1]][ind]
			acc = acc@y.values[[1]][ind]
			return(c(accuracy=acc, cutoff= cutoff))
		}))

		paucs = laply(r.preds, function(pred){
			pauc = performance(pred, "auc", fpr.stop= fpr.stop)
			pauc = pauc@y.values[[1]] / fpr.stop
			# print(pauc)
			return(pauc)
		}, .progress = "text")

		res[get.pred, ] = paucs

		##### running for smoothed data
		r.spreds = alply(spreds, 2, function(nd){
			pred <- prediction( nd, Y)			
		}, .progress = "text")

		saccs = t(laply(r.spreds, function(pred){
			acc = performance(pred, "acc")
			ind = which.max(acc@y.values[[1]])
			cutoff = acc@x.values[[1]][ind]
			acc = acc@y.values[[1]][ind]
			return(c(accuracy=acc, cutoff= cutoff))
		}))		

		spaucs = laply(r.spreds, function(pred){
			pauc = performance(pred, "auc", fpr.stop= fpr.stop)
			pauc = pauc@y.values[[1]] / fpr.stop
			# print(pauc)
			return(pauc)
		}, .progress = "text")

		sres[get.pred, ] = spaucs
		print(get.pred)
		predname = nii.stub(basename(fdf$img[get.pred]))
		predname = file.path(id.outdir, 
			paste0(predname, "_model_results", adder, ".Rda"))

		col.q = colQuantiles(preds, 
			probs = seq(0.95, 1, by=0.01))
		rownames(col.q) = colnames(preds)

	save(vol.data, vol.sdata, res, sres, lmod, benchmark,
		paucs, spaucs, fpr.stop,
		col.q,
		saccs, accs, brain.vol, not0100, 
		file = predname)
	print(correct)
}
# }
# }

# res = t(res)
# sres = t(sres)



