#################################################################
## This code is for plotting densities of ROI data
##
## Author: John Muschelli
## Last updated: May 20, 2014
#################################################################
#################################################################
rm(list=ls())
library(ggplot2)
library(fslr)
library(smallPDF)
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

correct = "none"
options = c("none", "N3", "N4", "N3_SS", "N4_SS",
		"SyN", "SyN_sinc", "Rigid", "Affine")

dval.df = hval.df = val.df = NULL


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



	##############################
	# Keeping files where predictors exist
	##############################
	outfiles = nii.stub(basename(fdf$img))
	outfiles = paste0(outfiles, "_predictors", adder, ".Rda")
	outfiles = file.path(fdf$outdir, outfiles)
	stopifnot(file.exists(outfiles))

	get.pred = 16
	x = fdf[get.pred,]

	##############################
	# Run lmod number of models - not all the models - leave out
	##############################
	mod.filename = file.path(outdir, 
		paste0("Collapsed_Models", adder, ".Rda"))
	load(mod.filename)
	nopred = seq(non.aggmods)
	
	fdf.run = fdf[-nopred,]
	nr = nrow(fdf)
	valid.ind = ceiling(nr/2)
	test.ind = seq( valid.ind +1, nr)
	valid.ind = seq(1, valid.ind)

	mygroup = "Validation"
	fdf.run = fdf.run[valid.ind, ]

	###################################
	# Create Data frames
	###################################
	for (get.pred in seq(nrow(fdf.run))){

		iddir = fdf.run$iddir[get.pred]
		id.outdir = fdf.run$outdir[get.pred]
		img.stub = nii.stub(fdf.run$img[get.pred], bn=TRUE)
		outname = img.stub
		outname = file.path(id.outdir, 
			paste0(outname, "_ROI_values", adder, ".Rda"))
		load(file=outname)
		dd = data.frame(x = dval$x, y = dval$y, img = img.stub,
			stringsAsFactors = FALSE)
		hd = data.frame(x = hval$mids, y = hval$density, 
			img = img.stub,
			stringsAsFactors = FALSE)
		vd = data.frame(x = vals, img = img.stub,
			stringsAsFactors = FALSE)	
		dval.df = rbind(dval.df, dd)
		hval.df = rbind(hval.df, hd)
		val.df = rbind(val.df, vd)
		print(get.pred)
	}


	###################################
	# Create ggplot objects
	###################################
	gd = ggplot(dval.df, aes(x=x, y=y, colour= img)) + 
		geom_line() + guides(colour=FALSE) + xlab("HU") +
		ggtitle(paste0("Type of Image Correction: ", correct, 
			", Group: ", mygroup))
	gh = gd %+% hval.df
	gv.all = ggplot(val.df, aes(x=x)) + 
		geom_line(stat = "density") + guides(colour=FALSE) + 
		xlab("HU") +
		ggtitle(paste0("Type of Image Correction: ", correct, 
			", Group: ", mygroup))
	gv = gv.all + aes(colour= img)

	cut.dval = dval.df[ dval.df$x <= 120, ]
	cut.hval = hval.df[ hval.df$x <= 120, ]
	cut.val = val.df[ val.df$x <= 120, ]
	gd.cut = gd %+% cut.dval
	gh.cut = gd %+% cut.hval
	gv.cut = gv %+% cut.val	
	gv.all.cut = gv.all %+% cut.val	

	###################################
	# Make pdf
	###################################
	pdfname = file.path(outdir, 
	    paste0("ROI_Density", adder, "_", mygroup, ".pdf"))
	# pdf(pdfname)
	if (file.exists(pdfname)) {
		file.remove(pdfname)
	}
    pdfobj = smallpdf()
		# print(gv)
		print(gv)
		print(gv.all)

		# print(gv)
		print(gv.cut)
		print(gv.all.cut)
    smallpdf.off(
    	pdfname = pdfname, 
        pdfobj = pdfobj, 
        clean = TRUE)
    # dev.off()
	print(correct)
}