####################################################################
## This code is for prediciton of Image Segmentation of CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
####################################################################
####################################################################
rm(list=ls())
library(plyr)
library(cttools)
library(ROCR)
library(matrixStats)
library(reshape2)
library(ggplot2)
library(fslr)
# library(car)
library(GGally)
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
correct = "SyN"
options = c("none", "N3", "N4", "N3_SS", "N4_SS",
		"SyN", "SyN_sinc", "Rigid", "Affine")
for (correct in options){
	
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
	

	short_predict = function(object, newdata, lthresh=  
		.Machine$double.eps^0.5){
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

	irow = 1
	x = fdf[irow,]


	outfiles = nii.stub(basename(fdf$img))
	outfiles = paste0(outfiles, "_model_results", adder, ".Rda")
	outfiles = file.path(fdf$outdir, outfiles)
	stopifnot(file.exists(outfiles))
	fdf = fdf[file.exists(outfiles), ]

	mod.filename = file.path(outdir, 
		paste0("Collapsed_Models", adder, ".Rda"))
	load(mod.filename)
	vol.sdatas = vol.datas = vol.data
	reses = sreses = res

	benches = vector(length= length(runpreds))
	all.saccs = all.accs = reses

	for (get.pred in 1:nrow(fdf)){

		idoutdir = fdf$outdir[get.pred]	
		predname = nii.stub(basename(fdf$img[get.pred]))
		predname = file.path(idoutdir, 
			paste0(predname, "_model_results", adder, ".Rda"))

		x = load(file = predname)
		reses[get.pred, ] = res[get.pred,]
		sreses[get.pred, ] = sres[get.pred,]

		all.accs[get.pred, ] = accs["accuracy",]
		all.saccs[get.pred, ] = saccs["accuracy",]

		
		vol.datas[get.pred, ] = vol.data[get.pred,]
		vol.sdatas[get.pred, ] = vol.sdata[get.pred,]

		benches[get.pred] = benchmark
		print(get.pred)
	}

	vol.data = vol.datas
	vol.sdata = vol.sdatas
	res = reses
	sres = sreses


	nopred = seq(non.aggmods)
	vd = vol.data
	vd = vd[-nopred,]
	nr = nrow(vd)
	valid.ind = ceiling(nr/2)
	test.ind = seq( valid.ind +1, nr)
	valid.ind = seq(1, valid.ind)

	outfile = file.path(outdir, 
		paste0("Model_performance_results", adder, ".Rda")
		)

	save(sres, res, vol.data, all.accs, all.saccs, benches,
		valid.ind, test.ind, nopred, lmod, 
		file = outfile)
	vsd = vol.sdata
	vsd = vsd[-nopred,]
	vsd = vsd[, !(colnames(vsd) %in% "truth")]

	###############################
	# Separate data into validation and test
	###############################
	truevol = vd[,"truth"]
	vd = vd[, !(colnames(vd) %in% "truth")]
	voldiff = vd - truevol
	adiff = abs(voldiff)

	sadiff = abs(vsd - truevol) 
	valid.svol = sadiff[valid.ind,]
	test.svol = sadiff[test.ind,]

	valid.vol = adiff[valid.ind,]
	test.vol = adiff[test.ind,]

	res = res[-nopred, ]
	sres = sres[-nopred, ]

	all.accs = all.accs[-nopred, ]
	all.saccs = all.saccs[-nopred, ]
	
	benches = benches[-nopred]

	valid.res = res[valid.ind,]
	test.res = res[test.ind,]

	valid.sres = sres[valid.ind,]
	test.sres = sres[test.ind,]

	valid.acc = all.accs[valid.ind,]
	valid.sacc = all.saccs[valid.ind,]

	test.acc = all.accs[test.ind,]
	test.sacc = all.saccs[test.ind,]

	valid.bench = benches[valid.ind]
	test.bench = benches[test.ind]

	bench.m = matrix(valid.bench, nrow=length(valid.bench), 
		ncol = ncol(valid.acc), byrow=FALSE)

	bench.diff = valid.acc - bench.m
	n_above = colSums(bench.diff > 0)

	above_bench = apply(bench.diff > 0, 2, all)
	best.acc = which(above_bench)
	biggest.acc_diff = which.max(colMeans(bench.diff))

	best.res= apply(valid.res, 1, which.max)
	best.sres= apply(valid.sres, 1, which.max)

	# sacc = valid.sacc
	# acc = valid.acc
	# inds = valid.ind
	# benches = valid.bench
	# vol = valid.vol
	# svol = valid.svol
	# sres = valid.sres
	ffdf = fdf[-nopred, ]
	group = "Test"
	if (group == "Training"){
		subset.ind = valid.ind
	}
	if (group == "Test"){
		subset.ind = test.ind
	}
	subset.ids = ffdf$id[subset.ind]


	pdfname = file.path(outdir, 
		paste0("Modeling_", group, "_Results", adder, ".pdf"))
	pdf(pdfname)
	#################################
	# Plot accuracy against the benchmarks
	#################################
	if ("mydf" %in% ls()){
		rm(list="mydf")
	}
	for (rundiff in c("Smoothed", "Unsmoothed")){
		if (rundiff  == "Smoothed") mydf = all.saccs
		if (rundiff  == "Unsmoothed") mydf = all.accs
		valids = data.frame(mydf)
		valids$id = ffdf$id
		# s = valid.sacc
		# colnames(s) = paste0("smooth_", colnames(s))
		# valids = cbind(valids, s)
		vbench = data.frame(benchmark=benches)
		vbench$id = ffdf$id
		long = melt(valids, id.vars = "id")
		long = merge(long, vbench, by="id", all=TRUE)

		long = long[ long$id %in% subset.ids, ]

		p = qplot(benchmark, value, colour = id, 
			data=long, facets = ~ variable) + guides(colour = FALSE) +
			geom_abline(intercept= 0, slope =1 ) + ylab("Accuracy") +
			xlab("Benchmark (predict all voxels 0)") + 
			ggtitle(rundiff)
		print(p)
		rm(list="mydf")
	}

	plotter = function(data, title="", ncol=4, nrow=5){
		long = melt(data)
		colnames(long) = c("id", "model", "value")
		probs = seq(0, 1, by=.1)
		quants = ddply(long, .(model), function(x) {
				quantile(x$value, probs = probs)
			})
		means = ddply(long, .(model), function(x) {
				c(mean=mean(x$value))
			})
		medians = ddply(long, .(model), function(x) {
				c(median=median(x$value))
			})	

		g = ggplot(data=long, aes(x=value, colour = model)) + 
			geom_density()
		g = ggplot(data=long, aes(x=value)) + geom_histogram() + 
			facet_wrap(~ model, nrow=nrow, ncol= ncol) 
		g = g + geom_vline(data=means, aes(xintercept=mean),
			colour="red") + 
			facet_wrap(~ model, nrow=nrow, ncol= ncol)
		g = g + geom_vline(data=medians, aes(xintercept=median),
			colour="green") + 
			facet_wrap(~ model, nrow=nrow, ncol= ncol)
		print(g + ggtitle(title))
		
		# par(mfrow = c(4, 4))
		# ddply(long, .(model), function(x) {
		# 	hist(x$value)
		# })

		return(quants)
	}



	volres = plotter(adiff[subset.ind, ], 
		title= "Difference in Predicted vs. True Volume, Unsmoothed")
	volsres = plotter(sadiff[subset.ind, ], 
		title= "Difference in Predicted vs. True Volume, Smoothed")


	add_params = function(g, type = "", params){
		plots = g$plots
		plots = unlist(plots)
		find = grep(paste0("ggally_", type), plots)	
		for (iplot in seq_along(find)){
			ifind = find[iplot]
			char = gsub(")$", ", ", g$plots[[ifind]])
			char = paste0(char, params, ")")
			g$plots[[ifind]] = char
		}
		return(g)
	}


	subtext = function(g, findtext = "", subber ="", ...){
		plots = g$plots
		plots = unlist(plots)
		for (iplot in seq_along(plots)){
			char = gsub(findtext, subber, g$plots[[iplot]], ...)
			g$plots[[iplot]] = char
		}
		return(g)
	}


	top5 = colMeans(sres[subset.ind,])
	top5 = names(sort(top5, decreasing=TRUE))[1:4]

	ggally_histogramDiag = function (data, mapping, ...) {
	    p <- ggplot(data, mapping) + 
	    scale_x_continuous() + 
	    scale_y_continuous() + 
	        geom_histogram(
	        	aes(y = ..density.. / max(..density..) * 
	        		diff(range(x, na.rm = TRUE))), ...)
	    p$type <- "diag"
	    p$subType <- "density"
	    p
	}


	ggally_abline = function (data, mapping, ...) 
	{
	    p <- ggplot(data = data, mapping)
	    p = p + geom_abline(...)
	    p <- p + geom_point(...)
	    p$type <- "continuous"
	    p$subType <- "smooth"
	    p
	}

	df = data.frame(sres[subset.ind, top5])
	g = ggpairs(df, 
		upper=list(continuous = "cor"), 
		lower = list(continuous = "smooth"),
		diag= list(continuous = "density")
		, axisLabels='show')

	# g = add_params(g, type="smooth", params = "col='blue'")
	# g2 = add_params(g, type="smooth", params = "colour='blue'")
	# g2 = subtext(g2, 
	# 	findtext = "ggally_densityDiag", 
	# 	subber = "ggally_histogramDiag")
	# g2

	# g3 = add_params(g, type="smooth", params = "colour='blue'")
	g3 = subtext(g, 
		findtext = "ggally_smooth", 
		subber = "ggally_abline")
	g3 = add_params(g3, type="abline", params = "intercept=0, slope=1")

	g3

	vres = plotter(res[subset.ind, top5], 
		title= "Partial AUC (under .1 FDR) Distribution, Unsmoothed", 
		ncol= 1, nrow=5)

	vsres = plotter(sres[subset.ind, top5], 
		title= "Partial AUC (under .1 FDR) Distribution, Smoothed",
		ncol= 1, nrow=5)

	both = sres[subset.ind, top5]
	colnames(both) = paste0("smooth_", colnames(both))
	both = cbind(both, res[subset.ind, top5])

	plotter(both, 
		title= "Partial AUC (under .1 FDR) Distribution",
		ncol= 1, nrow=8)

	dev.off()


}