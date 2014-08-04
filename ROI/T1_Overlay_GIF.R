rm(list=ls())
library(R.matlab)
library(oro.nifti)
library(plyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(knitrBootstrap)
library(xtable)
library(scales)
library(animation)
library(cttools)
library(fslr)
writepng = function(fname, tempimg= NULL, adder ='', rerun=TRUE, 
	text = "", col.y = alpha(hotmetal(10), 0.7), ...){
	pngname = gsub("\\.nii\\.gz", 
		paste0(adder, ".png"), fname)
	if (!file.exists(pngname) | rerun){
		pimg = readNIfTI(fname)
# mask.overlay		
		png(pngname, type="cairo")
		if (is.null(tempimg)) {
			orthographic(pimg, text=text, ...)
		} else {
			mask.overlay(tempimg, pimg, col.y=col.y, ...)
		}
		dev.off()
	}
	return(pngname)
}

img_cut = function(img, breaks, ...){
	cuts = cut(img, breaks=breaks, ...)
	# cuts = factor(cuts, levels)
	levs = levels(cuts)
	cuts = as.numeric(cuts)
	# res.p[ rs > ncut ] = cuts
	img@.Data = array(cuts, dim=dim(img))
	return(list(img=img, levs=levs))
}


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


#### baseline NIHSSS ata 
nihss = read.csv(file.path(basedir, "baseline_NIHSS.csv"), 
                 stringsAsFactors=FALSE)
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

outdir = file.path(basedir, "results")
whichdir = "reoriented"

# whichdir = "FLIRT"
rerun = FALSE

template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)

t.t1 = file.path(tempdir, "sct1_unsmooth.nii.gz")
temp.t1 = readNIfTI(t.t1)

setwd(progdir)


spm.binimg = file.path(outdir, paste0(whichdir, "_Binary_Sum_Image"))
spm.binimg = paste0(spm.binimg, ".nii.gz")


col.y = alpha(div_gradient_pal(low="blue", 
	mid=muted("red"), 
	high="yellow")(
	seq(0, 1, length=100)
	), 1)

fname = spm.binimg

spm2_t1_hot_overlay.png = writepng(spm.binimg, tempimg= temp.t1, 
	adder = "_t1_heat_overlay", window = c(500, 1000),
	col.y = col.y,
	rerun=TRUE, text ="Weighted Average")
tempimg = temp.t1
window = c(300, 1000)
text ="Population Average"
pimg = readNIfTI(fname)

#### 
breaks = seq(0, .45, by=.05)
col.cut = alpha(div_gradient_pal(low="blue", 
	mid="red", 
	high="yellow")(
	seq(0, 1, length=length(breaks)-1)
	), .7)

img = check_nifti(tempimg)
img.mask = check_nifti(pimg)
img@cal_min = window[1]
img@cal_max = window[2]
img[img < window[1]] = window[1]
img[img >= window[2]] = window[2]
img.mask[img.mask <= 0] = NA

clist = img_cut(pimg, breaks=breaks, include.lowest=FALSE)
cimg = clist$img
levs = clist$levs

plevs = levs
plevs = gsub("\\(", "", plevs)
plevs = gsub("\\]", "", plevs)
plevs = strsplit(plevs, ",")
plevs = lapply(plevs, as.numeric)
plevs = lapply(plevs, `*`, 100)
plevs = lapply(plevs, function(x){
	x[2] = x[2] - .01
	x
})
plevs = sapply(plevs, function(x){
	paste0(x[1], "-", x[2], "%")
})
fname = file.path(outdir, "Figure4_Percent.png")
png(fname, res=600, height=7, width=7, units="in")
ortho2(img, img.mask, col.y=col.cut, 
	ybreaks = breaks, 
	addlegend = TRUE,
	text="B", text.cex = 5,   
	text.x= 50, text.y=32,  
	leg.x = 5, leg.y= 60, 
	legend=plevs, 
	leg.col=col.cut, leg.cex=1.5,
	leg.title = "Percent of Sample\n with Hemorrhage")
dev.off()

fname = file.path(outdir, "Figure4_Proportion.png")
png(fname, res=600, height=7, width=7, units="in")
ortho2(img, img.mask, col.y=col.cut, 
	ybreaks = breaks, 
	addlegend = TRUE,
  text="B", text.cex = 5,   
	text.x= 50, text.y=32,  
	leg.x = 5, leg.y= 60, 
	legend=levs, 
	leg.col=col.cut, leg.cex=1.5,
	leg.title = "Proportion of Sample\n with Hemorrhage")
dev.off()

fname = file.path(outdir, "Figure4_Proportion_Final.png")
png(fname, res=600, height=7, width=7, units="in")
ortho2(img, img.mask, col.y=col.cut, 
       ybreaks = breaks, 
       addlegend = TRUE,
       leg.x = 14, leg.y= 60, 
       legend=levs, 
       leg.col=col.cut, leg.cex=1.5,
       leg.title = "Proportion of Sample\n with Hemorrhage")
dev.off()

# mask.overlay(tempimg, cimg, col.y=col.cut, text=text, window=window,
# 	ybreaks = seq(1, 11))





dtemp = dim(temp.t1)
xyz <- ceiling(dtemp/2)


setwd(outdir)

pimg = readNIfTI(spm.binimg)
adder = "_t1_heat_overlay"
vals = c(pimg@.Data)
df = data.frame(x = vals)
rm(list="vals")
df = df[ df$x > 0, , drop=FALSE]
g = ggplot(df, aes(x=x)) + geom_histogram(bin=.05) + 
	xlab("Proportion of patients with hemorrhage at voxel") +
	ylab("Number of Voxels") +  scale_y_continuous(labels = comma)	
g = g + theme(text= element_text(size=18)) + 
	ggtitle("Non-spatial Histogram of ICH Prevalence")	
g = g + annotate("text", 
					x = .4,
					y = 2e5,
					label = "A", size=20)


pdf(gsub("\\.nii\\.gz", "_histogram.pdf", spm.binimg))
	print(g)
dev.off()

png(gsub("\\.nii\\.gz", "_histogram.png", spm.binimg), 
    res=600, height=7, width=7, units="in")
	print(g) 
dev.off()


out.rda = gsub("\\.nii\\.gz", "_histogram_data.rda", spm.binimg)
save(df, file=out.rda)

pngname = gsub("\\.nii\\.gz", adder, spm.binimg)
movie.name = basename(paste0(pngname, ".gif"))
img.name = basename(pngname)



# myseq= range(which(pimg > 0, arr.ind=TRUE)[,3])
# myseq = seq(from=myseq[1], to=myseq[2])
# Z = dtemp[3]
# saveGIF({
# 	pb= txtProgressBar(min=0, max=length(myseq), style=3)
# 	for (iz in seq_along(myseq)){
# 		z = myseq[iz]
# 		adder = sprintf("_t1_heat_overlay_%03d", iz)
# 		pngname = gsub("\\.nii\\.gz", adder, fname)
# 		xyz[3] = z
# 		mask.overlay(tempimg, pimg, col.y=col.y, 
# 			window=window, text=text, xyz=xyz)
# 		setTxtProgressBar(pb, value=iz)

# 	}

# }, movie.name = movie.name, img.name= img.name)

# file.copy(file.path(tempdir(), movie.name), 
# 	file.path(outdir, movie.name), overwrite=TRUE)

# view.png(spm1_t1_hot_overlay.png)

