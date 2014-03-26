#################################
# Performs voxel-wise models for age/gender/volume
#################################
rm(list=ls())
library(oro.nifti)
library(plyr)
library(limma)
library(microbenchmark)
library(abind)
library(ggplot2)
library(data.table)
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
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

outdir = file.path(basedir, "results")

whichdir = "reoriented"

#### load voxel data
outfile = file.path(outdir, "Voxel_Matrix.Rda")
load(file=outfile )

#### keeping if over 10 people have ICH in that locaiton
ncut = 10
mat = mat[rs > ncut, ]

#### getting unique rows and duplications
dups = duplicated(mat)
dt = data.table(mat)
urows = unique(dt)

urows = t(urows)

# class(mat) = "numeric"
#### baseline NIHSSS ata 
nihss = read.csv(file.path(basedir, 
	"baseline_NIHSS.csv"), 
                 stringsAsFactors=FALSE)
nihss = nihss[ nihss$patientName %in% df$id, ]
nihss = nihss[ order(nihss$patientName), ]

nihss = nihss$nihss_total

keep = which(!is.na(nihss))
urows = urows[keep,]
nihss = nihss[keep]

rr = range(nihss)
scale.dens = function(x){
	d = density(x)
	## make width 1
	d$y = d$y / max(d$y) * 0.5
	d
}
nihss.list = alply(.data=urows, .margins = 2, function(x){
	c(g0 = scale.dens(nihss[!x]), g1=scale.dens(nihss[x]))
}, .progress= "text")


n = 50
l = nihss.list[[1]]
xlim = c(-0.5 + 0, 0.5 + n)
plot(y= l$g0.x, x= l$g0.y*.5, xlim=xlim, ylim = rr, type="l", 
	col="blue")
lines(y= l$g1.x, x= -l$g1.y* .5, col="red")

for (i in 2:n){
	l = nihss.list[[i]]

	lines(y= l$g0.x, x= l$g0.y*.5 + i-1, col="blue")
	lines(y= l$g1.x, x= -l$g1.y* .5+i-1, col="red")
}