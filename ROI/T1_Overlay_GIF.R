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
### need cairo for cluster
options(bitmapType='cairo')
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
source(file.path(progdir, "convert_DICOM.R"))
source(file.path(progdir, "fslhd.R"))
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


spm2_t1_hot_overlay.png = writepng(spm.binimg, tempimg= temp.t1, 
	adder = "_t1_heat_overlay", window = c(500, 1000),
	rerun=rr, text ="Weighted Average")
fname = spm.binimg
tempimg = temp.t1
window = c(500, 1000)
text ="Population Average"
col.y = alpha(hotmetal(10), 0.7)
dtemp = dim(temp.t1)
xyz <- ceiling(dtemp/2)


setwd(outdir)

adder = "_t1_heat_overlay"
pngname = gsub("\\.nii\\.gz", adder, fname)
movie.name = basename(paste0(pngname, ".gif"))
img.name = basename(pngname)
pimg = readNIfTI(fname)

vals = c(pimg@.Data)
df = data.frame(x = vals)
rm(list="vals")
df = df[ df$x > 0, , drop=FALSE]
g = ggplot(df, aes(x=x)) + geom_histogram() + 
	xlab("Proportion of patients with hemorrhage at voxel") +
	ylab("Frequency")
pdf(gsub("\\.nii\\.gz", ".pdf", fname))
	print(g)
dev.off()

myseq= range(which(pimg > 0, arr.ind=TRUE)[,3])
myseq = seq(from=myseq[1], to=myseq[2])
Z = dtemp[3]
saveGIF({
	pb= txtProgressBar(min=0, max=length(myseq), style=3)
	for (iz in seq_along(myseq)){
		z = myseq[iz]
		adder = sprintf("_t1_heat_overlay_%03d", iz)
		pngname = gsub("\\.nii\\.gz", adder, fname)
		xyz[3] = z
		mask.overlay(tempimg, pimg, col.y=col.y, 
			window=window, text=text, xyz=xyz)
		setTxtProgressBar(pb, value=iz)

	}

}, movie.name = movie.name, img.name= img.name)

file.copy(file.path(tempdir(), movie.name), 
	file.path(outdir, movie.name), overwrite=TRUE)

# view.png(spm1_t1_hot_overlay.png)

