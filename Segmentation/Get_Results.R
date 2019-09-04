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

segdir = file.path(progdir, "Segmentation")
source(file.path(segdir, "performance_functions.R"))

outdir = file.path(basedir, "results")
correct = "Rigid"
# options = c("none", "N3", "N4", "N3_SS", "N4_SS",
#         "SyN", "SyN_sinc", "Rigid", "Affine", "Rigid_sinc", 
#         "Affine_sinc")
# options = c("none", "N3_SS", "N4_SS", 
# 		"Rigid", "Rigid_sinc")
options = c("none", "Rigid")

# types = c("", "_include", "_zval", "_zval2")
# types = c("_zval2", "_zval_all", '_zval2_medztemp')
types = "_zval2"
# types = "_zval2"
# "_include_all", 
# types = "_include_all"
type = types[1]

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
		"Affine" = "_Affine",
		"Rigid_sinc" = "_Rigid_sinc",
		"Affine_sinc" = "_Affine_sinc")

	#### load voxel data
	outfile = file.path(outdir, "Voxel_Info.Rda")
	load(file=outfile )

	outfile = file.path(outdir, 
		"111_Filenames_with_volumes_stats.Rda")
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
		sens.cut.vol.datas = sens.cut.vol.sdatas = cut.vol.datas
		dice.cut.vol.datas = dice.cut.vol.sdatas = cut.vol.datas
		# pauc.cut.vol.tsdatas = pauc.cut.vol.sdatas

		benches = vector(length= length(runpreds))
		all.saccs = all.accs = reses
		all.dices = all.sdices = reses
		mod.saccs = mod.accs = reses
		mod.ssens = mod.sens = reses
		mod.sspec = mod.spec = reses
		mod.sdice = mod.dice = reses

		pauc.mod.saccs = pauc.mod.accs = reses
		pauc.mod.ssens = pauc.mod.sens = reses
		pauc.mod.sspec = pauc.mod.spec = reses
		pauc.mod.sdice = pauc.mod.dice = reses


		sens.mod.saccs = sens.mod.accs = reses
		sens.mod.ssens = sens.mod.sens = reses
		sens.mod.sspec = sens.mod.spec = reses
		sens.mod.sdice = sens.mod.dice = reses


		dice.mod.saccs = dice.mod.accs = reses
		dice.mod.ssens = dice.mod.sens = reses
		dice.mod.sspec = dice.mod.spec = reses
		dice.mod.sdice = dice.mod.dice = reses

		for (get.pred in 1:nrow(fdf)){
			
			# if (get.pred == 111) next;
			idoutdir = fdf$outdir[get.pred]	
			predname = nii.stub(basename(fdf$img[get.pred]))
			predname = file.path(idoutdir, 
				paste0(predname, "_model_results", adder, type,
					".Rda"))
			if (file.exists(predname)){
				x = load(file = predname)
				stopifnot(colnames(reses) == colnames(res))
				stopifnot(colnames(sreses) == colnames(sres))
				reses[get.pred, ] = res[get.pred,]
				sreses[get.pred, ] = sres[get.pred,]

				rownames(accs) = rownames(saccs) = 
					c("accuracy" ,"cutoff")

				all.accs[get.pred, ] = accs["accuracy",]
				all.saccs[get.pred, ] = saccs["accuracy",]

				all.dices[get.pred, ] = dices[1, ]
				all.sdices[get.pred, ] = sdices[1,]				

				vol.datas[get.pred, ] = vol.data[get.pred,]
				vol.sdatas[get.pred, ] = vol.sdata[get.pred,]

				cut.vol.datas[get.pred, ] = 
					cut.vol.data[get.pred,]
				cut.vol.sdatas[get.pred, ] = 
					cut.vol.sdata[get.pred,]
				cut.vol.tsdatas[get.pred, ] = 
					cut.vol.tsdata[get.pred,]

				sens.cut.vol.datas[get.pred, ] = 
					sens.cut.vol.data[get.pred,]
				sens.cut.vol.sdatas[get.pred, ] = 
					sens.cut.vol.sdata[get.pred,]

				dice.cut.vol.datas[get.pred, ] = 
					dice.cut.vol.data[get.pred,]
				dice.cut.vol.sdatas[get.pred, ] = 
					dice.cut.vol.sdata[get.pred,]


			####Acc cutoffs
			mod.accs[get.pred, ] = sapply(cut.tabs, get.acc)
			mod.saccs[get.pred, ] = sapply(cut.stabs, get.acc)

			mod.sens[get.pred, ] = sapply(cut.tabs, get.sens)
			mod.ssens[get.pred, ] = sapply(cut.stabs, get.sens)

			mod.spec[get.pred, ] = sapply(cut.tabs, get.spec)
			mod.sspec[get.pred, ] = sapply(cut.stabs, get.spec)	

			mod.dice[get.pred, ] = sapply(cut.tabs, get.dice)
			mod.sdice[get.pred, ] = sapply(cut.stabs, get.dice)	


			pauc.cut.vol.datas[get.pred, ] = 
				pauc.cut.vol.data[get.pred,]
			pauc.cut.vol.sdatas[get.pred, ] = 
				pauc.cut.vol.sdata[get.pred,]
			# pauc.cut.vol.tsdatas[get.pred, ] = 
			# 	pauc.cut.vol.tsdata[get.pred,]

			####pAUC cutoffs
			pauc.mod.accs[get.pred, ] = sapply(pauc.cut.tabs, 
				get.acc)
			pauc.mod.saccs[get.pred, ] = sapply(pauc.cut.stabs, 
				get.acc)

			pauc.mod.sens[get.pred, ] = sapply(pauc.cut.tabs, 
				get.sens)
			pauc.mod.ssens[get.pred, ] = sapply(pauc.cut.stabs, 
				get.sens)

			pauc.mod.spec[get.pred, ] = sapply(pauc.cut.tabs, 
				get.spec)
			pauc.mod.sspec[get.pred, ] = sapply(pauc.cut.stabs, 
				get.spec)	

			pauc.mod.dice[get.pred, ] = sapply(pauc.cut.tabs, 
				get.dice)
			pauc.mod.sdice[get.pred, ] = sapply(pauc.cut.stabs, 
				get.dice)

			####Dice cutoffs
			sens.mod.accs[get.pred, ] = sapply(sens.cut.tabs, 
				get.acc)
			sens.mod.saccs[get.pred, ] = sapply(sens.cut.stabs, 
				get.acc)

			sens.mod.sens[get.pred, ] = sapply(sens.cut.tabs, 
				get.sens)
			sens.mod.ssens[get.pred, ] = sapply(sens.cut.stabs,
				get.sens)

			sens.mod.spec[get.pred, ] = sapply(sens.cut.tabs, 
				get.spec)
			sens.mod.sspec[get.pred, ] = sapply(sens.cut.stabs,
				get.spec)	

			sens.mod.dice[get.pred, ] = sapply(sens.cut.tabs, 
				get.dice)
			sens.mod.sdice[get.pred, ] = sapply(sens.cut.stabs,
				get.dice)

			####Dice cutoffs
			dice.mod.accs[get.pred, ] = sapply(dice.cut.tabs, 
				get.acc)
			dice.mod.saccs[get.pred, ] = sapply(dice.cut.stabs,
				get.acc)

			dice.mod.sens[get.pred, ] = sapply(dice.cut.tabs, 
				get.sens)
			dice.mod.ssens[get.pred, ] = sapply(dice.cut.stabs,
				get.sens)

			dice.mod.spec[get.pred, ] = sapply(dice.cut.tabs, 
				get.spec)
			dice.mod.sspec[get.pred, ] = sapply(dice.cut.stabs,
				get.spec)	

			dice.mod.dice[get.pred, ] = sapply(dice.cut.tabs, 
				get.dice)
			dice.mod.sdice[get.pred, ] = sapply(dice.cut.stabs,
				get.dice)									

				benches[get.pred] = benchmark
				print(get.pred)
				rm(list=x[ x != "fpr.stop"])
			}
		}

		vol.data = vol.datas
		vol.sdata = vol.sdatas
		res = reses
		sres = sreses

		nopred = which(fdf$group == "Train")
		ffdf = fdf[-nopred, ]
		nr = nrow(ffdf)
		valid.ind = which(ffdf$group == "Validation")
		test.ind = which(ffdf$group == "Test")

		group = "Test"
		subset.ind = which(ffdf$group == group)

		varslice = ffdf$varslice[subset.ind]
		subset.ids = ffdf$id[subset.ind]
		all.truevol =  ffdf$truevol
		truevol = ffdf$truevol[subset.ind]
		full.subset.ind = match(subset.ids, fdf$id)

		outfile = file.path(outdir, 
			paste0("Model_performance_results", adder, type, 
				".Rda")
			)

		save(sres, res, valid.ind, test.ind,
			vol.datas, vol.sdatas,
			cut.vol.datas, cut.vol.sdatas,
			cut.vol.tsdatas,

			sens.cut.vol.datas, sens.cut.vol.sdatas,
			dice.cut.vol.datas, dice.cut.vol.sdatas,
			# pauc.cut.vol.tsdatas,
			pauc.cut.vol.datas, pauc.cut.vol.sdatas,
			# mod.tsaccs, mod.tssens, mod.tsspec,

			mod.saccs, mod.accs,
			mod.ssens, mod.sens,
			mod.sspec, mod.spec,
			mod.sdice, mod.dice,
	# sens cutoffs

			sens.mod.saccs, sens.mod.accs,
			sens.mod.ssens, sens.mod.sens,
			sens.mod.sspec, sens.mod.spec,
			sens.mod.sdice, sens.mod.dice,
	# dice cutoffs
			dice.mod.saccs, dice.mod.accs,
			dice.mod.ssens, dice.mod.sens,
			dice.mod.sspec, dice.mod.spec,
			dice.mod.sdice, dice.mod.dice,
	# pAUC cutoffs
			pauc.mod.saccs, pauc.mod.accs,
			pauc.mod.ssens, pauc.mod.sens,
			pauc.mod.sspec, pauc.mod.spec,
			pauc.mod.sdice, pauc.mod.dice,			

			all.accs, all.saccs, 
			all.dices, all.sdices, 

			benches,
			valid.ind, test.ind, nopred, 
			file = outfile)

		vsd = vol.sdata
		vsd = vsd[-nopred,]
		vsd = vsd[, !(colnames(vsd) %in% "truth")]

		###############################
		# Separate data into validation and test
		###############################
		vd = vol.data
		vd = vd[-nopred,]
		# truevol = vd[,"truth"]
		vd = vd[, !(colnames(vd) %in% "truth")]
		voldiff = vd - all.truevol
		adiff = abs(voldiff)


		######################################
		# Taking volume differences from truth and subsetting
		######################################
		####Accuracy cutoffs

		cut.diff = cut.vol.datas[-nopred, 
				!colnames(cut.vol.datas) %in% "truth"] - 
			all.truevol
		cut.adiff = abs(cut.diff)

		cut.sdiff = cut.vol.sdatas[-nopred, 
				!colnames(cut.vol.sdatas) %in% "truth"] -
			all.truevol
		cut.sadiff = abs(cut.sdiff)

		cut.tsdiff = cut.vol.tsdatas[-nopred, 
				!colnames(cut.vol.sdatas) %in% "truth"] -
			all.truevol
		cut.tsadiff = abs(cut.tsdiff)	

		####Sensitivity cutoffs
		sens.cut.diff = sens.cut.vol.datas[-nopred, 
				!colnames(sens.cut.vol.datas) %in% "truth"] - 
			all.truevol
		sens.cut.adiff = abs(sens.cut.diff)

		sens.cut.sdiff = sens.cut.vol.sdatas[-nopred, 
				!colnames(sens.cut.vol.sdatas) %in% "truth"] -
			all.truevol
		sens.cut.sadiff = abs(sens.cut.sdiff)

		####Dice cutoffs
		dice.cut.diff = dice.cut.vol.datas[-nopred, 
				!colnames(dice.cut.vol.datas) %in% "truth"] - 
			all.truevol
		dice.cut.adiff = abs(dice.cut.diff)

		dice.cut.sdiff = dice.cut.vol.sdatas[-nopred, 
				!colnames(dice.cut.vol.sdatas) %in% "truth"] -
			all.truevol
		dice.cut.sadiff = abs(dice.cut.sdiff)

		####pAUC cutoffs
		pauc.cut.diff = pauc.cut.vol.datas[-nopred, 
				!colnames(pauc.cut.vol.datas) %in% "truth"] - 
			all.truevol
		pauc.cut.adiff = abs(pauc.cut.diff)

		pauc.cut.sdiff = pauc.cut.vol.sdatas[-nopred, 
				!colnames(pauc.cut.vol.sdatas) %in% "truth"] -
			all.truevol
		pauc.cut.sadiff = abs(pauc.cut.sdiff)


		svoldiff = vsd - all.truevol
		sadiff = abs(svoldiff) 
		valid.svol = svoldiff[valid.ind,]
		test.svol = svoldiff[test.ind,]

		######################################
		# The ones with a's are absolute, 
		# should just do in plot code
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

		all.dices = all.dices[-nopred, ]
		all.sdices = all.sdices[-nopred, ]

		valid.acc = all.accs[valid.ind,]
		valid.sacc = all.saccs[valid.ind,]

		test.acc = all.accs[test.ind,]
		test.sacc = all.saccs[test.ind,]

		cut.all.accs = mod.accs[-nopred, ]
		cut.all.saccs = mod.saccs[-nopred, ]

		cut.all.dices = mod.dice[-nopred, ]
		cut.all.sdices = mod.sdice[-nopred, ]		

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

		cutoffs = c("", "sens.", "dice.", "pauc.")
		measures = c("accs", "sens", "dice", "spec")
		smooth = c("", "s")
		eg = expand.grid(cutoff = cutoffs, measure = measures,
			smooth = smooth, stringsAsFactors = FALSE)
		eg$obj = paste0(eg$cutoff, "mod.", eg$smooth, eg$measure)
		eg$best.value= NA
		eg$best = NA
		eg$best.median = eg$best.min = NA
		eg$best.worst = NA

		narm = TRUE
		for (iobj in seq(nrow(eg))){
			obj = get(eg$obj[iobj])
			obj = obj[-nopred, ]
			obj = obj[subset.ind, ]
			cn = colnames(obj)
			means = colMeans(obj, na.rm = narm)
			eg$best[iobj] = cn[which.max(means)]
			eg$best.value[iobj] = means[which.max(means)]

			meds = colMedians(obj, na.rm = narm)
			eg$best.median[iobj] = cn[which.max(meds)]

			mins = colMins(obj, na.rm = narm)
			eg$best.min[iobj] = cn[which.max(mins)]			
			eg$best.worst[iobj] = mins[which.max(mins)]
			# maxs = colMaxs(obj)
		}

		best.mat = ddply(eg, .(measure), function(x){
			best = which.max(x$best.value)
			x[best, c("obj", "best", "best.value")]
		})

		best.worst = ddply(eg, .(measure), function(x){
			best = which.max(x$best.worst)
			x[best, c("obj", "best", "best.worst")]
		})		

		vox.eg = eg
		######################################
		# Doing for volume differences
		cutoffs = c("", "sens.", "dice.", "pauc.")
		measures = c("diff", "adiff", 'novar.diff', 'novar.adiff')
		smooth = c("", "s")
		eg = expand.grid(cutoff = cutoffs, measure = measures,
			smooth = smooth, stringsAsFactors = FALSE)
		eg$obj = paste0(eg$cutoff, "cut.", eg$smooth, eg$measure)
		eg$best.value= NA
		eg$best = NA
		eg$best.median = eg$best.min = NA
		eg$best.worst = NA

		for (iobj in seq(nrow(eg))){
			novar = FALSE
			obj = eg$obj[iobj]
			ind = subset.ind
			if (grepl("novar", obj)){
				novar = TRUE
				ind = subset.ind[!varslice]
				obj = gsub("novar[.]", "", obj)
			}
			obj = get(obj)
			obj = obj[ind, ]
			cn = colnames(obj)
			means = colMeans(obj, na.rm = narm)
			eg$best[iobj] = cn[which.min(abs(means))]
			eg$best.value[iobj] = means[which.min(abs(means))]

			meds = colMedians(obj, na.rm = narm)
			eg$best.median[iobj] = cn[which.min(abs(meds))]

			maxs = colMaxs(abs(obj), na.rm = narm)
			eg$best.min[iobj] = cn[which.min(abs(maxs))]	
			eg$best.worst[iobj] = maxs[which.min(abs(maxs))]
			# maxs = colMaxs(obj)
		}

		vol.best.mat = ddply(eg, .(measure), function(x){
			best = which.min(abs(x$best.value))
			x[best, c("obj", "best", "best.value")]
		})

		vol.best.worst = ddply(eg, .(measure), function(x){
			best = which.min(abs(x$best.worst))
			x[best, c("obj", "best", "best.worst")]
		})

		vol.eg = eg
		best.mat = rbind(best.mat, vol.best.mat)
		best.worst = rbind(best.worst, vol.best.worst)

		outfile = file.path(outdir, 
			paste0("Best_", group, "_model_performance_results", 
				adder, type, 
				".Rda")
			)

		save(vox.eg, vol.eg, best.mat, best.worst, file=outfile)

		plotter = function(data, title="", varslice = NA,
			ncol=4, nrow=5,
			forcex=FALSE, xlimits=c(0, 1), 
			forcey = FALSE, ylimits = c(0, 13),
			tsize = 2, 
			tsize_all = 16,
			ylab="Frequency",
			loc.mult = 0.5,
			xlab = "value"){
			long = melt(data)
			colnames(long) = c("id", "model", "value")
			long$varslice = varslice
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
			mins = ddply(long, .(model), function(x) {
					c(ind.min=round(min(x$value), 3))
				})							
			maxes$max = means$max = medians$max = max(long$value)
			maxes$min = means$min = medians$min = min(long$value)

			mins$max = means$max 
			mins$min = means$min


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
			mins$height = height

			means$loc.mult = medians$loc.mult = 
				maxes$loc.mult = loc.mult
			mins$loc.mult = loc.mult

			g = ggplot(data=long, aes(x=value, 
				colour = model)) + 
				geom_density() 
			g2 = g + facet_wrap(~varslice)
			g = ggplot(data=long, aes(x=value)) + 
			geom_histogram() + 
				facet_wrap(~ model, nrow=nrow, ncol= ncol) 
			g2 = g + aes(fill=varslice)
			g = g + geom_vline(data=means, aes(xintercept=mean),
				colour="red") + 
				facet_wrap(~ model, nrow=nrow, ncol= ncol)
			g = g + geom_vline(data=medians, 
				aes(xintercept=median),
				colour="green") + 
				facet_wrap(~ model, nrow=nrow, ncol= ncol)
			g = g + geom_text(data=medians, 
				aes(x = (max + min) * loc.mult, 
				y=height, 
				label=paste0("Med: ", median)), size=tsize)
			g = g + geom_text(data=means, 
				aes(x = (max + min) * loc.mult, 
				y=height*.75, 
				label=paste0("Mean: ", mean)), size=tsize)
			g = g + geom_text(data=maxes, 
				aes(x = (max + min) * loc.mult, 
				y=height*.50, 
				label=paste0("Max: ", ind.max)), size=tsize)
			g = g + geom_text(data=mins, 
				aes(x = (max + min) * loc.mult, 
				y=height*.25, 
				label=paste0("Min: ", ind.min)), size=tsize)
			g = g + ggtitle(title) 
			g = g + ylab(ylab)
			g = g + xlab(xlab)
			g = g + theme(text = element_text(size = tsize_all))
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

		ssens = mod.ssens[subset.ind, ]
		sspec = mod.sspec[subset.ind, ]

		sens = mod.sens[subset.ind, ]
		spec = mod.spec[subset.ind, ]

		sens.ssens = sens.mod.ssens[subset.ind, ]
		sens.sspec = sens.mod.sspec[subset.ind, ]

		sens.sens = sens.mod.sens[subset.ind, ]
		sens.spec = sens.mod.spec[subset.ind, ]	

		melter = function(data, colname){
			x = melt(data)
			# colnames(x) = c("row.id", "model", colname)
			x[,3]
		}

		######################################
		### Sensitivity and Specificity of diff cutoffs
		######################################
		rows = melt(sens)
		rows = rows[,1:2]
		colnames(rows) = c("row.id", "model")
		allss = data.frame(sapply(list(ssens, sspec, sens, spec), 
			melter))
		colnames(allss) = c("acc.ssens", "acc.sspec", 
			"acc.sens", "acc.spec")
		ss = data.frame(sapply(list(sens.ssens, sens.sspec, 
			sens.sens, sens.spec), 
			melter))
		colnames(ss) = c("sens.ssens", "sens.sspec", 
			"sens.sens", "sens.spec")

		allss = cbind(rows, allss, ss)
		long = melt(allss, id.vars = c("row.id", "model"))
		long$type = gsub("(.*)[.].*", "\\1", long$variable)
		long$variable = gsub("(.*)[.](.*)", "\\2", long$variable)
		long$smooth = grepl("^ss", long$variable)
		long$variable = gsub("^ss", "s", long$variable)
		rm(list=c("sens", "spec", "ssens", "sspec", "ss",
			"sens.ssens", "sens.sspec", "sens.sens", "sens.spec"))
		wide = reshape(data=long, direction = "wide", 
			idvar = c("row.id", "model", "type", "smooth"), 
			timevar= "variable")
		colnames(wide) = gsub("value[.]", "", colnames(wide))

		ord = wide[order(-wide$spec, wide$sens), ]
		ord = ord[ ord$model %in% c("mod_agg"), ]

		g = ggplot(ord, aes(x=1-spec, y=sens, colour=model)) +
			geom_line() +
			facet_wrap( ~ type + smooth) 

		run = wide[ ! wide$model %in% paste0("model", 1:10), ]

		stats = ddply(run, .(model, type, smooth), summarise,
			mean_sens = mean(sens, na.rm=TRUE),
			min_sens = min(sens, na.rm=TRUE),
			mean_spec = mean(spec, na.rm=TRUE),
			min_spec = min(spec, na.rm=TRUE)
			)


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
		# L = list()
		for (icut in c("cut", "")){
			for (rundiff in c("Smoothed", "Unsmoothed")){
				if (icut =="") {
					if (rundiff  == "Smoothed") {
						mydf = all.saccs
					}
					if (rundiff  == "Unsmoothed") {
						mydf = all.accs
					}
				}
				if (icut =="cut") {
					if (rundiff  == "Smoothed") {
						mydf = cut.all.saccs
					}
					if (rundiff  == "Unsmoothed") {
						mydf = cut.all.accs
					}
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
					c("mod_agg", "min", "gam", "rf"), ])
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

		# ba.plotter(voldiff[subset.ind, ], truth = truevol,
#   title= "Difference in Predicted vs. True Volume, Unsmoothed")

		# ba.plotter(svoldiff[subset.ind, ], 
#   title= "Difference in Predicted vs. True Volume, Smoothed")

		###########################################
		# Volume differences and BA plots for best models
		###########################################
		plotter(voldiff[subset.ind, top5], 
			varslice =fdf$varslice[subset.ind],
		title=
		"Difference in Predicted vs. True Volume, Unsmoothed")

		plotter(svoldiff[subset.ind, stop5], 
		title=
			"Difference in Predicted vs. True Volume, Smoothed")

		plotter(cut.sdiff[subset.ind, stop5], 
		title=
		"Difference in Predicted (Cut) vs. True Volume, Smoothed")


		ba.plotter(voldiff[subset.ind, top5], truth = truevol,
		title=
			"Difference in Predicted vs. True Volume, Unsmoothed")

		ba.plotter(svoldiff[subset.ind, stop5], truth = truevol,
		title=
			"Difference in Predicted vs. True Volume, Smoothed")


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
	# "Difference in Predicted vs. True Volume, Unsmoothed, Cut")

	# volres = plotter(cut.adiff[subset.ind, 
	# 	c("mod_agg", "median", "gmean", "gam", "min")], 
	#   title= 
	# "Difference in Predicted vs. True Volume, Unsmoothed, Cut")
	# volsres = plotter(cut.sadiff[subset.ind, ], 
	#   title= 
	# "Difference in Predicted vs. True Volume, Smoothed, Cut")	

	# volsres = plotter(cut.sadiff[subset.ind, 
	# 	c("mod_agg", "median", "gmean", "gam", "min")], 
	#   title= 
	#   "Difference in Predicted vs. True Volume, Smoothed, Cut")

	# volsres = plotter(sadiff[subset.ind, top5], 
	#title= "Difference in Predicted vs. True Volume, Smoothed")

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
				char = gsub(findtext, subber, 
					g$plots[[iplot]], ...)
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
			forcey = TRUE, ylimits = c(0, 20))
		dev.off()
		

		message("After modeling_auc\n")
		mydf = cut.all.saccs
		##################
		# Model aggregate plots
		###################
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

		long = long[ long$variable %in% "mod_agg", ]

		p = qplot(benchmark, value, colour = id, 
			data=long) + 
			guides(colour = FALSE) +
			geom_abline(intercept= 0, slope =1 ) + 
			ylab("Accuracy") +
			xlab("Benchmark (predict all voxels 0)") + 
			ggtitle(paste0(rundiff, " ", icut))
		
		pngname = file.path(outdir, 
			paste0("Modeling_", group, "_AccPlot", adder, 
				type, ".png"))
		png(pngname, type= "cairo",  res=600,  width = 7, 
			height = 7, units = "in")		
		print(p)
		dev.off()
		rm(list="mydf")



		run.vd = cut.sdiff[subset.ind, 'mod_agg', drop=FALSE]
		colnames(run.vd) = "Smoothed Aggregate Model"
		pngname = file.path(outdir, 
			paste0("Modeling_", group, "_VolDiff", adder, 
				type, "_Final.png"))	
		png(pngname, type= "cairo",  res=600,  width = 7, 
			height = 7, units = "in")					
		plotter(run.vd, 
			title= "Difference in Predicted vs. True Volume")
		dev.off()


		pickmod = "mod_agg"

		for (pickmod in c("rf", "mod_agg", "gam")){
			stubber = function(inname){
				paste0("Modeling_", group, "_", inname, 
					adder, type, "_Final", pickmod)
			}

		
			both = cbind("Aggregate Model" = res[subset.ind, 
				pickmod],
				"Smoothed Aggregate Model" = sres[subset.ind, 
				pickmod])

			pngname = file.path(outdir, 
				paste0(stubber("AUC"), ".png")
				)
			png(pngname, type= "cairo",  res=600,  width = 7, 
				height = 7, units = "in")
				plotter(both, 
				title= paste0("Partial AUC (under ", fpr.stop, 
				" FDR) Distribution"),
				ncol= 1, nrow=2, forcex = TRUE, 
				loc.mult = 0.1,
				tsize = 5, xlab = "pAUC" )
			dev.off()


			both = cbind("Aggregate Model" = 
				dice.mod.dice[full.subset.ind, pickmod],
				"Aggregate Model after Smoothing Probabilities" = 
				dice.mod.sdice[full.subset.ind, pickmod])

			pngname = file.path(outdir, 
				paste0(stubber("Dice"), ".png")
				)
			png(pngname, type= "cairo",  res=600,  width = 7, 
				height = 7, units = "in")
				plotter(both, 
				title= paste0("Dice Similarity Index Distribution ", 
					"\n in ", group ," Scans"),
				ncol= 1, nrow=2, forcex = TRUE, 
				loc.mult = 0.1,
				tsize_all = 16,
				tsize = 5, xlab = "Dice Similarity Index" )
			dev.off()

			both = cbind("Aggregate Model" = 
				mod.dice[full.subset.ind, pickmod],
				"Aggregate Model after Smoothing Probabilities" = 
				mod.sdice[full.subset.ind, pickmod])

			pngname = file.path(outdir, 
				paste0(stubber("Dice"), 
					"_Other_Crit.png")
				)
			png(pngname, type= "cairo",  res=600,  width = 7, 
				height = 7, units = "in")
				plotter(both, 
				title= paste0("Dice Similarity Index Distribution ", 
					"\n in ", group ," Scans"),
				ncol= 1, nrow=2, forcex = TRUE, 
				loc.mult = 0.1,
				tsize_all = 16,
				tsize = 5, 
				xlab = "Dice Similarity Index" )
			dev.off()		

			both = cbind(
				"Smoothed Aggregate Model Predictions" = 
				dice.mod.sdice[full.subset.ind, pickmod])

			pngname = file.path(outdir, 
				paste0(stubber("Dice"), 
					"_Smooth.png")
				)
			png(pngname, type= "cairo",  res=600,  width = 7, 
				height = 7, units = "in")
				plotter(both, 
				title= paste0("Dice Similarity Index Distribution ", 
					"\n in ", group ," Scans"),
				ncol= 1, nrow=2, forcex = TRUE, 
				loc.mult = 0.1,
				tsize_all = 16,
				tsize = 5, 
				xlab = "Dice Similarity Index" )
			dev.off()		


			df = data.frame(pred = vsd[subset.ind, pickmod])
			df$true = truevol
			df$avg = (df$true + df$pred)/2
			df$diff = df$pred - df$true
			df$varslice = varslice
			pngname = file.path(outdir, 
				paste0(stubber("BA_Volume_Plot"), 
					"_Smooth.png")
				)			
			png(pngname, type= "cairo",  res=600,  width = 7, 
				height = 7, units = "in")
			qplot(x=avg, 
				y=diff, 
				data=df, 
				geom=c("point")) + 
				ggtitle("Bland-Altman Plot of ICH Volume") + 
				xlab("Average (Predicted + Manual)/2 ICH Volume") + 
				geom_smooth(se=FALSE) +
				theme(text = element_text(size = 16)) + 
				ylab("Difference (Predicted - Manual) ICH Volume") + 
				annotate("text", 
					label = "Manual > Predicted ICH Volume", x = 30, 
					y = -15, size = 5, colour = "black") +
				annotate("text", 
					label = "Predicted > Manual ICH Volume", x = 30, 
					y = 20, size = 5, colour = "black")
			dev.off()
		} # end pickmod

		print(type)
	} # type loop
} # correct loop