#################################
# Gets outcomes for voxelwise regressions
#################################
rm(list=ls())
library(oro.nifti)
library(plyr)
library(limma)
library(microbenchmark)
library(scales)
library(RColorBrewer)
library(data.table)
library(cttools)
library(fslr)
library(smallpdf)
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


#### load voxel data
outfile = file.path(outdir, "Voxel_Matrix.Rda")
load(file=outfile )

template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)

t1 = file.path(tempdir, "betsct1_unsmooth.nii.gz")
temp.t1 = readNIfTI(t1)

dtemp = dim(temp)

stopifnot(prod(dtemp) == nrow(mat))


outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile)


thresholds = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1)

tt = temp
tt[ !is.na(tt) ] = NA
tt[ is.na(tt) ] = TRUE
inds = which(tt > 0, arr.ind=TRUE)
inds = inds[rs > ncut, ]



for (outcome in c("GCS", "NIHSS")){
	whichdir = "reoriented"
	adder = paste0(outcome, "_")
	if (outcome == "NIHSS"){
		adder = ""
	}
	print(outcome)


	outfile = file.path(outdir, paste0(adder, "Pvalue_Matrix.Rda"))
	load(file=outfile)

	#### keeping if over 10 people have ICH in that locaiton


	device = "png"

	overlayimg = temp.t1
	addto = ""
	if (all(overlayimg == temp.t1)){
		addto = "_t1"
	}

	runmods = c(1, 2, 5, 6)
	addtexts = LETTERS[1:4]
	iimod = 1
	for (iimod in seq_along(runmods)){

		imod = runmods[iimod]

		fp = file.path(outdir, 
			paste0(adder, "Regression_Map_heatcol", imod, 
				addto, "_Final.", device))
		open_dev(fp, res=600, height=7, width=7, units = "in")
			res.p = temp
			res.p[!is.na(res.p)] = NA
			y = res[,"X",imod]
			stopifnot(all(y > 0))
			cuts = cut(y, breaks=thresholds, include.lowest=FALSE)
			cuts = factor(cuts, levels=rev(levels(cuts)))
			levs = levels(cuts)
			cuts = as.numeric(cuts)
			res.p[ rs > ncut ] = cuts
			# res.p[rs > ncut] = -log(y)
			# res.p[rs > ncut] = as.numeric(y <= 0.05)
			# res.p[ rs > ncut ]  = 1-y
			cols = c("blue", "deepskyblue1", "green", "yellow", 
				"orange", "red")

			##### all this depends on what unique cuts are shown in the
			#### perspective
			xyz = c(91, 109, 91)

			cuts1 = cuts[ which(inds[,1] == xyz[1])]
			cuts2 = cuts[ which(inds[,2] == xyz[2])]
			cuts3 = cuts[ which(inds[,3] == xyz[3])]
			ucuts = sort(unique(c(cuts1, cuts2, cuts3)))
			# cols = cols[ucuts]
			# ucuts = c(0, 6)
			# cols = cols[ucuts]		
			cols = alpha(cols, .7)
			addtext = addtexts[iimod]
			mask.overlay(overlayimg, res.p, window = c(0, 1000),
				col.y=cols, zlim.y = c(0, 6),
				addlegend = TRUE,
				leg.x = 15, leg.y= 45, 
				legend=rev(levs), 
				leg.col=alpha(rev(cols), 1), leg.cex=1.5,
				leg.title = "P-value", text=addtext, text.cex =5,
				text.x= 32, text.y=58
			)
			# col.y = alpha(hotmetal(10), 0.7)	
			# mask.overlay(temp, res.p, col.y=col.y)
			# col.y[1] = alpha("white", 0.7)
			# mask.overlay(temp, res.p, col.y=col.y)
			# cols = alpha(brewer.pal(9, "Reds"), 0.7)
			# mask.overlay(temp, res.p, col.y=cols)
			# rm(list="res.p")
			
		dev.off()
		# view.png(fp)
	# }

	}


}