###############################
# make matrix of p-values for the group
###############################
rm(list=ls())
library(oro.nifti)
library(plyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(knitrBootstrap)
library(xtable)
library(scales)
library(fslr)
library(cttools)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
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
dtemp = dim(temp)


maskfile = file.path(tempdir, "scct_mask.nii.gz")
tempmask = readNIfTI(maskfile)




#### get all roi images
imgs = list.files(pattern='ROI.*.nii', 
	full.names=TRUE, path=basedir, recursive=TRUE)

df = data.frame(dir = dirname(imgs), fname = basename(imgs),
	stringsAsFactors = FALSE)
df = df[ !grepl("2mm_", df$fname), ]
if (grepl("reoriented", whichdir)) {
	df = df[ grepl("bws", df$fname), ]
}
df = df[ !grepl("affine9", df$fname), ]

df$id = df$dir
df$id = gsub(basedir, "", df$id)
df$id = gsub(".*(\\d\\d\\d-(\\d|)\\d\\d\\d)(.*)", "\\1", df$id)
df$id = as.numeric(gsub("-", "", df$id))

### grab only the first image
df = ddply(df, .(id), function(x){
	x[1,]
})

nimgs = nrow(df)
dn = vector("list", length=2)
dn[[1]] = 1:prod(dtemp)
dn[[2]] = df$fname
read.img = function(fname){
	img = readNIfTI(fname)
	img[is.nan(img) | is.na(img)] = 0
	# res = which(img > 0)
	# return(res)
	return(c(img > 0))
}
iimg = 1

df$img = gsub("bws", "w", df$fname)
df$img = gsub("ROI\\.", ".", df$img)

n = 21

img = readNIfTI(file.path(df$dir[1], df$img[1]))
roi = readNIfTI(file.path(df$dir[1], df$fname[1]))

masked = img
masked[tempmask == 0] = 0
masked[masked > 100 | masked < 0] = 0

diffimg = img
diffimg@.Data = (img - temp)/temp
diffimg = cal_img(diffimg)

diff2 = img
diff2@.Data = (img - temp)
diff2 = cal_img(diff2)

mdiff = diffimg
mdiff[tempmask == 0] = 0
mdiff = cal_img(mdiff)

z_z = zscore_img(img, mask=tempmask)

dat = data.frame(int=c(masked), temp=c(temp),
	rat=c(mdiff), roi=c(roi), diff= c(diff2), z=c(z_z))
dat = dat[ c(tempmask > 0), ]
dat = dat[ dat$int > 20 & dat$int < 90, ]
diffcut = 5
dat = dat[ dat$rat < diffcut, ]
dat = dat[ dat$rat > 0, ]
dat = dat[ dat$diff >= 10, ]
dat = dat[ dat$z > 0, ]


pred = mdiff
pred[ img < 20 | img > 90 ] = 0
pred[ mdiff > diffcut] = 0
pred[ mdiff < 0] = 0
pred[ diff2 < 10] = 0
pred[ z_z < 0] = 0
pred = cal_img(pred)

samp.ind = sample(nrow(dat), size=1e5)
ggplot(dat[samp.ind,], aes(x=int, y=diff, colour=factor(roi))) + 
	geom_point()

ggplot(dat[samp.ind,], aes(x=diff, y=z, colour=factor(roi))) + 
	geom_point()
# ggplot(dat[samp.ind,], aes(x=temp, y=diff, colour=roi)) + geom_point()

