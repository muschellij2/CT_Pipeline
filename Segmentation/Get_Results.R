##################################################################
## This code is for prediciton of Image Segmentation of CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
##################################################################
##################################################################
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
segdir = file.path(progdir, "Segmentation")
basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

outdir = file.path(basedir, "results")
correct = "none"
options = c("none", "N3", "N4", "N3_SS", "N4_SS",
		"SyN", "SyN_sinc", "Rigid", "Affine")

types = c("", "_include", "_zval", "_zval2")
# "_include_all", 
# types = "_include_all"
type = types[4]


for (correct in options){

	print(correct)
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
	    X <- model.matrix(Terms, m, 
	    	contrasts.arg = object$contrasts)
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


	for (type in types){
		outfiles = nii.stub(basename(fdf$img))
		outfiles = paste0(outfiles, "_model_results", adder, type,
			".Rda")
		outfiles = file.path(fdf$outdir, outfiles)
		# stopifnot(file.exists(outfiles))
		# fdf = fdf[file.exists(outfiles), ]

		mod.filename = file.path(outdir, 
			paste0("Collapsed_Models", adder, ".Rda"))
		load(mod.filename)
		vol.sdatas = vol.datas = vol.data
		reses = sreses = res

		cut.vol.datas = cut.vol.sdatas = vol.data
		cut.vol.tsdatas = cut.vol.sdatas

		pauc.cut.vol.datas = pauc.cut.vol.sdatas = cut.vol.datas
		# pauc.cut.vol.tsdatas = pauc.cut.vol.sdatas

		benches = vector(length= length(runpreds))
		all.saccs = all.accs = reses
		mod.saccs = mod.accs = reses
		mod.ssens = mod.sens = reses
		mod.sspec = mod.spec = reses

		get.acc = function(tab){
			sum(diag(tab)) / sum(tab)
		}
		get.sens = function(tab){
			tab["TRUE", "TRUE"] / sum(tab[, "TRUE"])
		}
		get.spec = function(tab){
			tab["FALSE", "FALSE"] / sum(tab[, "FALSE"])
		}

		for (get.pred in 1:nrow(fdf)){
			
			if (get.pred == 111) next;
			idoutdir = fdf$outdir[get.pred]	
			predname = nii.stub(basename(fdf$img[get.pred]))
			predname = file.path(idoutdir, 
				paste0(predname, "_model_results", adder, type,
					".Rda"))

			x = load(file = predname)
			reses[get.pred, ] = res[get.pred,]
			sreses[get.pred, ] = sres[get.pred,]

			all.accs[get.pred, ] = accs["accuracy",]
			all.saccs[get.pred, ] = saccs["accuracy",]

			vol.datas[get.pred, ] = vol.data[get.pred,]
			vol.sdatas[get.pred, ] = vol.sdata[get.pred,]

			cut.vol.datas[get.pred, ] = cut.vol.data[get.pred,]
			cut.vol.sdatas[get.pred, ] = cut.vol.sdata[get.pred,]
			cut.vol.tsdatas[get.pred, ] = cut.vol.tsdata[get.pred,]

			mod.accs[get.pred, ] = sapply(cut.tabs, get.acc)
			mod.saccs[get.pred, ] = sapply(cut.stabs, get.acc)
			# mod.tsaccs[get.pred, ] = sapply(cut.tabs.smooth, get.acc)

			mod.sens[get.pred, ] = sapply(cut.tabs, get.sens)
			mod.ssens[get.pred, ] = sapply(cut.stabs, get.sens)
			# mod.tssens[get.pred, ] = sapply(cut.tabs.smooth, get.sens)

			mod.spec[get.pred, ] = sapply(cut.tabs, get.spec)
			mod.sspec[get.pred, ] = sapply(cut.stabs, get.spec)		
			# mod.tsspec[get.pred, ] = sapply(cut.tabs.smooth, get.spec)

			pauc.cut.vol.datas[get.pred, ] = 
				pauc.cut.vol.data[get.pred,]
			pauc.cut.vol.sdatas[get.pred, ] = 
				pauc.cut.vol.sdata[get.pred,]
			# pauc.cut.vol.tsdatas[get.pred, ] = 
			# 	pauc.cut.vol.tsdata[get.pred,]
					

			benches[get.pred] = benchmark
			print(get.pred)
			rm(list=x[ x != "fpr.stop"])
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
			paste0("Model_performance_results", adder, type, 
				".Rda")
			)

		save(sres, res, valid.ind, test.ind,
			vol.data, vol.sdata,
			cut.vol.datas, cut.vol.sdatas,
			cut.vol.tsdatas,
			# pauc.cut.vol.tsdatas,
			pauc.cut.vol.datas, pauc.cut.vol.sdatas,
			mod.saccs, mod.accs,
			# mod.tsaccs, mod.tssens, mod.tsspec,
			mod.sens, mod.ssens,
			mod.spec, mod.sspec,
			all.accs, all.saccs, benches,
			valid.ind, test.ind, nopred, 
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


		######################################
		# Taking volume differences from truth and subsetting
		######################################
		cut.diff = cut.vol.datas[-nopred, 
				!colnames(cut.vol.datas) %in% "truth"] - 
			truevol
		cut.adiff = abs(cut.diff)

		cut.sdiff = cut.vol.sdatas[-nopred, 
				!colnames(cut.vol.sdatas) %in% "truth"] -
			truevol
		cut.sadiff = abs(cut.sdiff)

		cut.tsdiff = cut.vol.tsdatas[-nopred, 
				!colnames(cut.vol.sdatas) %in% "truth"] -
			truevol
		cut.tsadiff = abs(cut.tsdiff)				

		svoldiff = vsd - truevol
		sadiff = abs(svoldiff) 
		valid.svol = svoldiff[valid.ind,]
		test.svol = svoldiff[test.ind,]

		######################################
		# The ones with a's are absolute, should just do in plot code
		######################################
		valid.vol = voldiff[valid.ind,]
		test.vol = voldiff[valid.ind,]

		valid.avol = adiff[valid.ind,]
		test.avol = adiff[test.ind,]

		valid.savol = sadiff[valid.ind,]
		test.savol = sadiff[test.ind,]

		cut.valid.svol = cut.sdiff[valid.ind,]
		cut.test.svol = cut.sdiff[test.ind,]

		cut.valid.savol = cut.sadiff[valid.ind,]
		cut.test.savol = cut.sadiff[test.ind,]

		cut.valid.tsvol = cut.tsdiff[valid.ind,]
		cut.test.tsvol = cut.tsdiff[test.ind,]	

		cut.valid.tsavol = cut.tsadiff[valid.ind,]
		cut.test.tsavol = cut.tsadiff[test.ind,]

		cut.valid.vol = cut.diff[valid.ind,]
		cut.test.vol = cut.diff[test.ind,]		

		cut.valid.avol = cut.adiff[valid.ind,]
		cut.test.avol = cut.adiff[test.ind,]

		######################################
		# Getting pAUC
		######################################
		res = res[-nopred, ]
		sres = sres[-nopred, ]

		valid.res = res[valid.ind,]
		test.res = res[test.ind,]

		valid.sres = sres[valid.ind,]
		test.sres = sres[test.ind,]

		######################################
		# Getting Accuracy
		######################################
		all.accs = all.accs[-nopred, ]
		all.saccs = all.saccs[-nopred, ]

		valid.acc = all.accs[valid.ind,]
		valid.sacc = all.saccs[valid.ind,]

		test.acc = all.accs[test.ind,]
		test.sacc = all.saccs[test.ind,]

		cut.all.accs = mod.accs[-nopred, ]
		cut.all.saccs = mod.saccs[-nopred, ]

		cut.valid.acc = cut.all.accs[valid.ind,]
		cut.valid.sacc = cut.all.saccs[valid.ind,]

		cut.test.acc = mod.accs[test.ind,]
		cut.test.sacc = mod.saccs[test.ind,]

		######################################
		# Getting Benchmarks
		######################################
		benches = benches[-nopred]
		valid.bench = benches[valid.ind]
		test.bench = benches[test.ind]	

		bench.m = matrix(valid.bench, nrow=length(valid.bench), 
			ncol = ncol(valid.acc), byrow=FALSE)

		######################################
		# Getting Difference in accuracy and benchmark
		######################################
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
		group = "Training"
		if (group == "Training"){
			subset.ind = valid.ind
		}
		if (group == "Test"){
			subset.ind = test.ind
		}
		subset.ids = ffdf$id[subset.ind]
		truth = truevol[subset.ind]


		plotter = function(data, title="", 
			ncol=4, nrow=5,
			forcex=FALSE, xlimits=c(0, 1), 
			forcey = FALSE, ylimits = c(0, 13)){
			long = melt(data)
			colnames(long) = c("id", "model", "value")
			probs = seq(0, 1, by=.1)
			quants = ddply(long, .(model), function(x) {
					quantile(x$value, probs = probs)
				})
			means = ddply(long, .(model), function(x) {
					c(mean=round(mean(x$value), 3))
				})
			medians = ddply(long, .(model), function(x) {
					c(median=round(median(x$value), 3))
				})	
			maxes = ddply(long, .(model), function(x) {
					c(ind.max=round(max(x$value), 3))
				})			
			maxes$max = means$max = medians$max = max(long$value)
			maxes$min = means$min = medians$min = min(long$value)


			fx = function(x) {
				rx = range(x)
				bw = diff(rx)/30
				bins = seq(from=rx[1], to=rx[2], by=bw)
				h = hist(x, breaks=bins, plot=FALSE)
				mc = max(h$counts)
			}
			height=max(daply(long, .(model), function(x){
				fx(x$value)
				# print(x)
			}))
			means$height = height
			medians$height = height
			maxes$height = height


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
			g = g + geom_text(data=medians, aes(x = (max + min)/2, 
				y=height, 
				label=paste0("Med: ", median)), size=2)
			g = g + geom_text(data=means, aes(x = (max + min)/2, 
				y=height*.8, 
				label=paste0("Mean: ", mean)), size=2)
			g = g + geom_text(data=maxes, aes(x = (max + min)/2, 
				y=height*.6, 
				label=paste0("Max: ", ind.max)), size=2)
			g = g + ggtitle(title)
			if (forcex) g= g + xlim(xlimits)
			if (forcey) g= g + ylim(ylimits)
			print(g)
			

			# par(mfrow = c(4, 4))
			# ddply(long, .(model), function(x) {
			# 	hist(x$value)
			# })

			return(quants)
		}

		ba.plotter = function(data, truth, 
			title="", 
			ncol=4, nrow=5,
			forcex=FALSE, xlimits=c(0, 1), 
			forcey = FALSE, ylimits = c(0, 13)){
			long = melt(data)
			colnames(long) = c("id", "model", "value")
			long$truth = truth

			g = ggplot(data=long, aes(x=truth, y=value)) + 
			geom_point() + geom_smooth(se=FALSE) +
				facet_wrap(~ model, nrow=nrow, ncol= ncol) 
			g = g + ggtitle(title)
			if (forcex) g= g + xlim(xlimits)
			if (forcey) g= g + ylim(ylimits)
			print(g)

			return(invisible())
		}

		pdfname = file.path(outdir, 
			paste0("Modeling_", group, "_Results", adder, 
				type, ".pdf"))
		pdf(pdfname)
		#################################
		# Plot accuracy against the benchmarks
		#################################
		if ("mydf" %in% ls()){
			rm(list="mydf")
		}
		for (icut in c("cut", "")){
			for (rundiff in c("Smoothed", "Unsmoothed")){
				if (icut =="") {
					if (rundiff  == "Smoothed") mydf = all.saccs
					if (rundiff  == "Unsmoothed") mydf = all.accs
				}
				if (icut =="cut") {
					if (rundiff  == "Smoothed") mydf = cut.all.saccs
					if (rundiff  == "Unsmoothed") mydf = cut.all.accs
				}			

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
					data=long, facets = ~ variable) + 
					guides(colour = FALSE) +
					geom_abline(intercept= 0, slope =1 ) + 
					ylab("Accuracy") +
					xlab("Benchmark (predict all voxels 0)") + 
					ggtitle(paste0(rundiff, " ", icut))
				print(p)
				print(p %+% long[ long$variable %in% 
					c("mod_agg", "min", "gmean", "gam"), ])
				rm(list="mydf")
			}
		}


		nkeep = 4

		stop5 = colMeans(sres[subset.ind,])
		stop5 = names(sort(stop5, decreasing=TRUE))[1:nkeep]

		top5 = colMeans(res[subset.ind,])
		top5 = names(sort(top5, decreasing=TRUE))[1:nkeep]


		###########################################
		# Volume differences and BA plots for all models
		###########################################
		# plotter(voldiff[subset.ind, ], 
		#   title= "Difference in Predicted vs. True Volume, Unsmoothed")

		# plotter(svoldiff[subset.ind, ], 
		#   title= "Difference in Predicted vs. True Volume, Smoothed")

		# ba.plotter(voldiff[subset.ind, ], truth = truth,
		#   title= "Difference in Predicted vs. True Volume, Unsmoothed")

		# ba.plotter(svoldiff[subset.ind, ], 
		#   title= "Difference in Predicted vs. True Volume, Smoothed")

		###########################################
		# Volume differences and BA plots for best models
		###########################################
		plotter(voldiff[subset.ind, top5], 
		  title= "Difference in Predicted vs. True Volume, Unsmoothed")

		plotter(svoldiff[subset.ind, stop5], 
		  title= "Difference in Predicted vs. True Volume, Smoothed")

		ba.plotter(voldiff[subset.ind, top5], truth = truth,
		  title= "Difference in Predicted vs. True Volume, Unsmoothed")

		ba.plotter(svoldiff[subset.ind, stop5], truth = truth,
		  title= "Difference in Predicted vs. True Volume, Smoothed")


		###########################################
		# Absolute Volume differences for all models
		###########################################
		# volres = plotter(adiff[subset.ind, ], 
		#   title= paste0("Abs Difference in Predicted vs. ", 
		#   	"True Volume, Unsmoothed"))	
		# volsres = plotter(sadiff[subset.ind, ], 
		#   title= paste0("Abs Difference in Predicted vs. ", 
		#   	"True Volume, Smoothed"))

		###########################################
		# Absolute Volume differences for best models
		###########################################
		volres = plotter(adiff[subset.ind, top5], 
		  title= paste0("Abs Difference in Predicted vs. ", 
		  	"True Volume, Unsmoothed"),
		  ncol= 1, nrow=nkeep)	
		volsres = plotter(sadiff[subset.ind, stop5], 
		  title= paste0("Abs Difference in Predicted vs. ", 
		  	"True Volume, Smoothed"),
		  ncol= 1, nrow=nkeep)	


		# volres = plotter(cut.adiff[subset.ind, ], 
		#   title= 
		#   	"Difference in Predicted vs. True Volume, Unsmoothed, Cut")

		# volres = plotter(cut.adiff[subset.ind, 
		# 	c("mod_agg", "median", "gmean", "gam", "min")], 
		#   title= 
		#   	"Difference in Predicted vs. True Volume, Unsmoothed, Cut")	
		# volsres = plotter(cut.sadiff[subset.ind, ], 
		#   title= 
		#   "Difference in Predicted vs. True Volume, Smoothed, Cut")	

		# volsres = plotter(cut.sadiff[subset.ind, 
		# 	c("mod_agg", "median", "gmean", "gam", "min")], 
		#   title= 
		#   "Difference in Predicted vs. True Volume, Smoothed, Cut")		

		# volsres = plotter(sadiff[subset.ind, top5], 
		# 	title= "Difference in Predicted vs. True Volume, Smoothed")

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


		ggally_abline = function (data, mapping, ...) {
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
		g3 = add_params(g3, type="abline", 
			params = "intercept=0, slope=1")

		g3

		###########################################
		# pAUC for best models
		###########################################
		plotter(res[subset.ind, top5], title= 
			paste0("Partial AUC (under ", fpr.stop, 
			" FDR) Distribution, Unsmoothed"), 
			ncol= 1, nrow=nkeep, forcex = TRUE)

		plotter(sres[subset.ind, stop5], title= 
			paste0("Partial AUC (under ", fpr.stop, 
			" FDR) Distribution, Smoothed"),
			ncol= 1, nrow=nkeep, forcex = TRUE)

		both = sres[subset.ind, stop5]
		colnames(both) = paste0("smooth_", colnames(both))
		both = cbind(both, res[subset.ind, top5])

		plotter(both, title= paste0("Partial AUC (under ", 
			fpr.stop, 
			" FDR) Distribution"),
			ncol= 1, nrow=nkeep*2, forcex = TRUE)

		dev.off()


		pngname = file.path(outdir, 
			paste0("Modeling_", group, "_AUC", adder, 
				type, ".png"))
		png(pngname, type= "cairo",  res=600,  width = 7, 
			height = 7, units = "in")
			plotter(both, 
			title= paste0("Partial AUC (under ", fpr.stop, 
			" FDR) Distribution"),
			ncol= 1, nrow=nkeep*2, forcex = TRUE,
			forcey = TRUE, ylimits = c(0, 13))
		dev.off()
		print(type)
	} # type loop
} # correct loop