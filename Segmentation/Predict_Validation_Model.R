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
library(getopt)
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

segdir = file.path(progdir, "Segmentation")
source(file.path(segdir, "performance_functions.R"))

outdir = file.path(basedir, "results")

correct = "none"
# options = c("none", "N3", "N4", "N3_SS", "N4_SS",
#         "SyN", "SyN_sinc", "Rigid", "Affine", "Rigid_sinc", 
#         "Affine_sinc")
options = c("none", "N3_SS", "N4_SS",
	"Rigid",  "Rigid_sinc")
# options = c("none", "Rigid_sinc")
# options = "Rigid"

spec = matrix(c(
	'correct', 'c', 1, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

correct = opt$correct
print(opt)

my.tab <- function(
  x, 
  y, 
  dnames=c("x", "y")) {
  x = as.numeric(x)
  y = as.numeric(y)
  stopifnot(all(unique(c(x,y)) %in% c(0, 1, NA)))
  tt = sum(x * y)
  t1=sum(x)
  t2=sum(y)
  tab = matrix(c(length(x)-t1-t2+tt,  t1-tt, t2-tt, tt), 2, 2)
  n = list(c("FALSE", "TRUE"), c("FALSE", "TRUE"))
  names(n) = dnames
  dimnames(tab) = n
  tab = as.table(tab)
  return(tab) 
}

types = c("_zval2", "_zval_all", "_zval2_medztemp")
# , "_zval2"ls
# "_include_all", 
type = types[1]


keep.obj = ls()

# for (correct in options){
	
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
		"Affine" = "_Affine",
		"Rigid_sinc" = "_Rigid_sinc",
		"Affine_sinc" = "_Affine_sinc")

	#### load voxel data
	outfile = file.path(outdir, "Voxel_Info.Rda")
	load(file=outfile )

    outfile = file.path(outdir, 
    	"111_Filenames_with_volumes_stats.Rda")
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

	cut.vol.data = cut.vol.sdata = cut.vol.tsdata = vol.data
	pauc.cut.vol.data = pauc.cut.vol.sdata = cut.vol.data
	pauc.cut.vol.tsdata = pauc.cut.vol.sdata
	sens.cut.vol.data = sens.cut.vol.sdata = cut.vol.data
	dice.cut.vol.data = dice.cut.vol.sdata = cut.vol.data

    # train = all.df[samps,]
    # test = all.df[!samps,]

    fname = file.path(outdir, 
        paste0("Validation_data", adder, ".Rda"))
    x = load(fname)

    all.df$mode = fdf$mode[match(mult.df$img, fdf$img)]

	stopifnot(all(all.df$Y %in% c(0, 1)))
	
	set.seed(20150225)
	all.df$subset = FALSE
	all.df$subset[sample(nrow(all.df), size= 1e6)] = TRUE

	df = all.df[ all.df$subset, ]
	df = cbind(df, mult.df[ all.df$subset, ])
	######################################
	# Get volume of ROI
	######################################	
	Y =  df$Y

	fpr.stop = .01

	cc = complete.cases(df)
	stopifnot(all(cc))

	preds = t(laply(all.mods, short_predict, newdata= df, 
              .progress = "text"))

	#### we only think values 0 to 100 are actually blood
	rowMean = rowMeans(preds)
	rowMed = rowMedians(preds)
	rowMax = rowMaxs(preds)
	rowMin = rowMins(preds)
	rowProd = exp(rowSums(log(preds)))
    nr = nrow(preds)
    rowGeom = rowProd ^ (1/nr)

	colnames(preds) = colnames(res)[seq(ncol(preds))]


	test.gam.pred = rep(0, length=nrow(df))
	

	gam.pred = predict(gam.mod, 
		df, 
		newdata.guaranteed = TRUE,
		type="response", block.size=1e5)
    cat("GAM Prediction \n")
    
    gam.pred = as.numeric(gam.pred)
    gam.pred[gam.pred > 1] = 1
	test.gam.pred[ df$subset ] = gam.pred


	# for (imod in  seq(lmod)){
	# 	irow = which(vol.diff$imod == imod & 
	# 		vol.diff$get.pred == get.pred)
	preds = cbind(preds, rowMean, rowMed, rowMax, rowMin,
		rowProd, rowGeom, gam=test.gam.pred)
	rownames(preds) = NULL

	#### if the mask did not contain these voxels, 
	#### then they are 0

	preds[ !df$multiplier, ] = 0

	cut.filename = file.path(outdir, 
	paste0("Model_Cutoffs", adder, ".Rda"))

	x = load(file=cut.filename)

	scut.filename = file.path(outdir, 
	paste0("Smooth_Model_Cutoffs", adder, ".Rda"))

	y = load(file=scut.filename)

	stopifnot(all(names(all.cuts) == colnames(preds)))
	stopifnot(all(names(all.pauc.cuts) == colnames(preds)))

	stopifnot(all(names(all.scuts) == colnames(preds)))
	stopifnot(all(names(all.spauc.cuts) == colnames(preds)))

	xpreds = preds

	for (type in types){
		if (type == ""){
			df$multiplier = 1
		}
		if (type == "_include"){
			df$multiplier = df$include
		}	
		if (type == "_include_all"){
			df$multiplier = df$include.all
		}
		if (type == "_zval"){
			df$multiplier = df$zval
		}		
		if (type == "_zval2"){
			df$multiplier = df$zval2
		}	
		if (type == "_zval_all"){
			df$multiplier = df$zval_all
		}				
		if (type == "_zval2_medztemp"){
		  df$multiplier = df$zval2_medztemp
		}				    
		preds = xpreds
		#### only values 0 to 100 are likely to be ROI
		# preds = preds * df$in0100
		preds = preds * df$multiplier


		################################
		# Running with variable slice thickness
		################################

		run_vol = function(runpred){
			img = nim
			img[is.na(img)] = 0
			img@.Data[keep.ind] = runpred
			vol = get_roi_vol(img, dcmtables)$truevol 
		}

		rpreds = t(t(preds) > all.cuts)
		# cut.vol.data[get.pred, ] = c(vol.roi, vols)


		rpreds = t(t(preds) > all.pauc.cuts)
		# pauc.cut.vol.data[get.pred, ] = c(vol.roi, vols)

		rpreds = preds	
		# vol.data[get.pred, ] = c(vol.roi, vols)


		rpreds = t(t(preds) > all.sens.cuts)
		# sens.cut.vol.data[get.pred, ] = c(vol.roi, vols)

		rpreds = t(t(preds) > all.dice.cuts)
		# dice.cut.vol.data[get.pred, ] = c(vol.roi, vols)		

		# pred.imgs = vector(mode="list", length=ncol(preds))
		# spreds = matrix(nrow=nrow(preds), ncol=ncol(preds))
		# colnames(spreds) = colnames(preds)
		# pb = txtProgressBar(max=ncol(preds), style=3)
		# for (ipred in seq(ncol(preds))){
		# 	nd = preds[, ipred]
		# 	img = nim
		# 	nd[roi.not.in] = 0
		# 	img@.Data[keep.ind] = nd
		# 	img = cal_img(img)
		# 	img[is.na(img)]= 0
		# 	img = datatyper(img, 
		# 		datatype= convert.datatype()$FLOAT32, 
		# 		bitpix = convert.bitpix()$FLOAT32)
		# 	cn = colnames(preds)[ipred]
		# 	# pred.imgs[[ipred]] = img
		# 	outimg = nii.stub(fdf$img[get.pred], bn=TRUE)
		# 	outimg = file.path(id.outdir, 
		# 		paste0(outimg, "_", cn, adder))	
		# 	writeNIfTI(img, filename = outimg )

		# 	#### smooth the results
		# 	sm.img  = mean_image(img, nvoxels = 1)
		# 	# sm.img = local_moment(img, nvoxels = 1, moment =1,
		# 	# 	center = FALSE)[[1]]
		# 	sm.img[abs(sm.img) <  .Machine$double.eps^0.5 ] = 0
		# 	sm.img = sm.img[keep.ind]
		# 	spreds[, ipred] = sm.img
		# 	spreds[roi.not.in, ipred] = 0

		# 	sm.img = nim
		# 	sm.img@.Data[keep.ind] = spreds[, ipred]
		# 	sm.img[is.na(sm.img)]= 0
		# 	sm.img[roi.not.in] = 0

		# 	sm.img = cal_img(sm.img)
		# 	sm.img  = datatyper(sm.img , 
		# 		datatype= convert.datatype()$FLOAT32, 
		# 		bitpix = convert.bitpix()$FLOAT32)

		# 	outimg = paste0(outimg, "_smoothed", adder, type)
		# 	writeNIfTI(sm.img, filename = outimg )
		# 	setTxtProgressBar(pb, ipred)

		# }
		# close(pb)

		##################
		# Deleted for space issues
		# predname = nii.stub(basename(fdf$img[get.pred]))
		# predname = file.path(outdir, 
		# 	paste0(predname, "_predictions", adder, ".Rda"))
		# save(preds, spreds, Y, file=predname, compress=TRUE)
		##################



		# rpreds = t(t(spreds) > all.scuts)
		# if (correct %in% c("none", "N3_SS", "N4_SS")) {
		# 	vols = aaply(rpreds, 2, run_vol,
		# 		.progress = "text")
		# } else {
		# 	vols = colSums(rpreds)
		# 	vols = vols * vres
		# }
		# cut.vol.sdata[get.pred, ] = c(vol.roi, vols)

		# rpreds = t(t(spreds) > all.cuts)
		# if (correct %in% c("none", "N3_SS", "N4_SS")) {
		# 	vols = aaply(rpreds, 2, run_vol,
		# 		.progress = "text")
		# } else {
		# 	vols = colSums(rpreds)
		# 	vols = vols * vres
		# }
		# cut.vol.tsdata[get.pred, ] = c(vol.roi, vols)

		# rpreds = t(t(spreds) > all.spauc.cuts)
		# if (correct %in% c("none", "N3_SS", "N4_SS")) {
		# 	vols = aaply(rpreds, 2, run_vol,
		# 		.progress = "text")
		# } else {
		# 	vols = colSums(rpreds)
		# 	vols = vols * vres
		# }
		# pauc.cut.vol.sdata[get.pred, ] = c(vol.roi, vols)

		# rpreds = t(t(spreds) > all.pauc.cuts)
		# if (correct %in% c("none", "N3_SS", "N4_SS")) {
		# 	vols = aaply(rpreds, 2, run_vol,
		# 		.progress = "text")
		# } else {
		# 	vols = colSums(rpreds)
		# 	vols = vols * vres
		# }
		# pauc.cut.vol.tsdata[get.pred, ] = c(vol.roi, vols)

		# rpreds = spreds
		# if (correct %in% c("none", "N3_SS", "N4_SS")) {
		# 	vols = aaply(rpreds, 2, run_vol,
		# 		.progress = "text")
		# } else {
		# 	vols = colSums(rpreds)
		# 	vols = vols * vres
		# }			
		# vol.sdata[get.pred, ] = c(vol.roi, vols)


		# rpreds = t(t(spreds) > all.ssens.cuts)
		# if (correct %in% c("none", "N3_SS", "N4_SS")) {
		# 	vols = aaply(rpreds, 2, run_vol,
		# 		.progress = "text")
		# } else {
		# 	vols = colSums(rpreds)
		# 	vols = vols * vres
		# }
		# sens.cut.vol.sdata[get.pred, ] = c(vol.roi, vols)		

		# rpreds = t(t(spreds) > all.sdice.cuts)
		# if (correct %in% c("none", "N3_SS", "N4_SS")) {
		# 	vols = aaply(rpreds, 2, run_vol,
		# 		.progress = "text")
		# } else {
		# 	vols = colSums(rpreds)
		# 	vols = vols * vres
		# }
		# dice.cut.vol.sdata[get.pred, ] = c(vol.roi, vols)

		runtabs = function(x, y){
			tabs = vector(mode="list", length=ncol(x))
			for (itab in seq(ncol(x))){
				tabs[[itab]] = my.tab(x[, itab] > y[itab], df$Y)
				print(itab)
			}
			names(tabs) = colnames(x)
			tabs
		}
		cut.tabs = runtabs(preds, all.cuts)
		# cut.stabs = runtabs(spreds, all.scuts)
		# cut.tabs.smooth = runtabs(spreds, all.cuts)

		pauc.cut.tabs = runtabs(preds, all.pauc.cuts)
		# pauc.cut.stabs = runtabs(spreds, all.spauc.cuts)
		# pauc.cut.tabs.smooth = runtabs(spreds, all.pauc.cuts)

		sens.cut.tabs = runtabs(preds, all.sens.cuts)
		# sens.cut.stabs = runtabs(spreds, all.ssens.cuts)
		# sens.cut.tabs.smooth = runtabs(spreds, all.sens.cuts)	

		dice.cut.tabs = runtabs(preds, all.dice.cuts)
		# dice.cut.stabs = runtabs(spreds, all.sdice.cuts)
		# dice.cut.tabs.smooth = runtabs(spreds, all.dice.cuts)

		
		run_acc = function(nd){
			pred <- prediction( nd, Y)	
			acc = get_acc(pred) 

			pauc = get_pauc(pred, fpr.stop)
			dice = get_max_dice(pred)
			list(acc = acc, pauc = pauc, dice = dice)
		}

		r.res = alply(preds, 2, run_acc, .progress = "text")

		accs = sapply(r.res, `[[`, "acc")
		paucs = sapply(r.res, `[[`, "pauc")
		dices = sapply(r.res, `[[`, "dice")
		names(paucs) = colnames(accs) =
			colnames(dices) =
			as.character(unlist(attr(r.res,"split_labels")))

		##### running for smoothed data
		# s.res = alply(spreds, 2, run_acc, .progress = "text")
		# saccs = sapply(s.res, `[[`, "acc")
		# spaucs = sapply(s.res, `[[`, "pauc")


	# 	col.q = colQuantiles(preds, 
	# 		probs = seq(0.95, 1, by=0.01))
	# 	rownames(col.q) = colnames(preds)

	# 	save(vol.data, vol.sdata, 
	# 		cut.vol.data, cut.vol.sdata,
	# 		pauc.cut.vol.data, pauc.cut.vol.sdata,
	# 		pauc.cut.vol.tsdata,
	# 		sens.cut.vol.data, sens.cut.vol.sdata,
	# 		dice.cut.vol.data, dice.cut.vol.sdata,
	# 		cut.vol.tsdata, 
	# 		cut.tabs.smooth,
	# 		pauc.cut.tabs.smooth,
	# 		res, sres, lmod, benchmark,
	# 		paucs, spaucs, fpr.stop,
	# 		pauc.cut.tabs, all.cuts, all.pauc.cuts,
	# 		pauc.cut.stabs, 
	# 		sens.cut.tabs, sens.cut.stabs, sens.cut.tabs.smooth,
	# 		dice.cut.tabs, dice.cut.stabs, dice.cut.tabs.smooth,
	# 		cut.tabs, 
	# 		cut.stabs, 	
	# 		col.q,
	# 		saccs, accs, brain.vol, not0100, 
	# 		file = predname)
	# 	rm(list=c("s.res", "r.res", "spreds"))
	# }
	# 	print(correct)
# }
# }
# }

# res = t(res)
# sres = t(sres)



