#####################################
## Author: John Muschelli
## Date: January 20, 2014
## Purpose: Read in the AAL atlas labels and make R objects
## that can be used later for overlap metrics.
#####################################
#####################################
rm(list=ls())
library(R.matlab)
library(oro.nifti)
library(plyr)
library(XML)
library(stringr)
homedir = "/Applications"
basedir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  basedir = "/dexter/disk2/smart/stroke_ct/ident"
}
spmdir = file.path(homedir, "spm8")

tempdir = file.path(basedir, "Template")
outdir = file.path(tempdir, "atlases")

fsldir <- system("echo $FSLDIR", intern=TRUE)
if (fsldir == "") {
  cmd <- paste("FSLDIR=/usr/local/fsl; FSLOUTPUTTYPE=NIFTI_GZ; ", 
    "export FSLDIR FSLOUTPUTTYPE; sh ${FSLDIR}/etc/fslconf/fsl.sh;")
}

atlas_dir = file.path(fsldir, "data", "atlases")

### this is to just truncate the last row and column because the 
### CT Template is this dimensions
outdim = c(181, 217, 181)

cutdown = function(img, outdim= c(181, 217, 181)){
	x = img@.Data 
	x = x[1:outdim[1], 1:outdim[2], 1:outdim[3]]
	img@.Data = x
	# dim(img) = outdim
	img@dim_ = c(3, outdim, 1, 1, 1, 1)
	return(img)
}

#### extracting Talairach Image an labels###########

atlas = "Talairach"
tatlas_dir = file.path(atlas_dir, atlas)
xmlfile = file.path(atlas_dir, paste0(atlas, ".xml"))

xx = xmlParse(xmlfile)
indices = xpathSApply(xx, "/atlas/data/label", xmlGetAttr, "index")
labels = xpathSApply(xx, "/atlas/data/label", xmlValue)
labs = strsplit(labels, "\\.")
labs = lapply(labs, str_trim)
labs = do.call("rbind", labs)

df = data.frame(labs, stringsAsFactors=FALSE)
colnames(df) = c("Area", "Lobe", "Lobule", "Tissue_Type", "Broadmann")
df = rbind(rep("Uncategorized", ncol(df)), df)
df$index = c(0, as.numeric(indices) + 1)

tal.df = df


img = readNIfTI(
	file.path(tatlas_dir, 
		"Talairach-labels-1mm.nii.gz"))
uimg = sort(unique(c(img)))
all.ind = c(0, df$index)
stopifnot(all(uimg %in% all.ind))

img = cutdown(img)

tal.img = img


#### extracting MNI Image an labels###########


atlas = "MNI"
tatlas_dir = file.path(atlas_dir, atlas)
xmlfile = file.path(atlas_dir, paste0(atlas, ".xml"))

xx = xmlParse(xmlfile)
indices = xpathSApply(xx, "/atlas/data/label", xmlGetAttr, "index")
labels = xpathSApply(xx, "/atlas/data/label", xmlValue)
labs = str_trim(labels)

df = data.frame(labs, stringsAsFactors=FALSE)
colnames(df) = "Label"
df = rbind(rep("Uncategorized", ncol(df)), df)
df$index = c(0, as.numeric(indices) + 1)

mni.df = df


img = readNIfTI(
	file.path(tatlas_dir, 
		"MNI-maxprob-thr0-1mm.nii.gz"))
uimg = sort(unique(c(img)))
all.ind = c(0, df$index)
stopifnot(all(uimg %in% all.ind))

img = cutdown(img)

mni.img = img


#### extracting Harvard-Oxford Cortical Image an labels###########

atlas = "HarvardOxford"
tatlas_dir = file.path(atlas_dir, atlas)
xmlfile = file.path(atlas_dir, paste0(atlas, "-Cortical.xml"))

xx = xmlParse(xmlfile)
indices = xpathSApply(xx, "/atlas/data/label", xmlGetAttr, "index")
labels = xpathSApply(xx, "/atlas/data/label", xmlValue)
labs = strsplit(labels, ",")
labs = lapply(labs, function(x) {
	x = str_trim(x)
	x = c(x, rep(NA, 2-length(x)))
	return(x)
})
labs = do.call("rbind", labs)


df = data.frame(labs, stringsAsFactors=FALSE)
colnames(df) = c("Label", "Division")
df = rbind(rep("Uncategorized", ncol(df)), df)
df$index = c(0, as.numeric(indices) + 1)

hoxcort.df = df

img = readNIfTI(
	file.path(tatlas_dir, 
		"HarvardOxford-cort-maxprob-thr0-1mm.nii.gz"))
uimg = sort(unique(c(img)))
all.ind = c(0, df$index)
stopifnot(all(uimg %in% all.ind))

img = cutdown(img)

hoxcort.img = img


#### extracting Harvard-Oxford Subcortical Image an labels###########

atlas = "HarvardOxford"
tatlas_dir = file.path(atlas_dir, atlas)
xmlfile = file.path(atlas_dir, paste0(atlas, "-Subcortical.xml"))

xx = xmlParse(xmlfile)
indices = xpathSApply(xx, "/atlas/data/label", xmlGetAttr, "index")
labels = xpathSApply(xx, "/atlas/data/label", xmlValue)
labs = str_trim(labels)


df = data.frame(labs, stringsAsFactors=FALSE)
colnames(df) = c("Label")
df = rbind(rep("Uncategorized", ncol(df)), df)
df$index = c(0, as.numeric(indices) + 1)

hoxsubcort.df = df
hoxsubcort.df$Label = gsub("Ventrical", "Ventricle", 
	hoxsubcort.df$Label)

img = readNIfTI(
	file.path(tatlas_dir, 
		"HarvardOxford-sub-maxprob-thr0-1mm.nii.gz"))
uimg = sort(unique(c(img)))
all.ind = c(0, df$index)
stopifnot(all(uimg %in% all.ind))

img = cutdown(img)

hoxsubcort.img = img

#### extracting JHU Eve atlas Type I an labels###########

atlas = "JHU_MNI_SS_WMPM_Type-I"
tatlas_dir = file.path(tempdir)
txtfile = file.path(tempdir, paste0(atlas, "_SlicerLUT.txt"))

xx = read.table(txtfile)
xx = xx[, 1:2]
colnames(xx) = c("index", "Label")
jhut1.df = xx

img = readNIfTI(
	file.path(tatlas_dir, 
		paste0(atlas, ".nii.gz")))
uimg = sort(unique(c(img)))
all.ind = jhut1.df$index
stopifnot(all(uimg %in% all.ind))

jhut1.img = img




#### extracting JHU Eve atlas Type II an labels###########

atlas = "JHU_MNI_SS_WMPM_Type-II"
tatlas_dir = file.path(tempdir)
txtfile = file.path(tempdir, paste0(atlas, "_SlicerLUT.txt"))

xx = read.table(txtfile)
xx = xx[, 1:2]
colnames(xx) = c("index", "Label")
jhut2.df = xx

img = readNIfTI(
	file.path(tatlas_dir, 
		paste0(atlas, ".nii.gz")))
uimg = sort(unique(c(img)))
all.ind = jhut2.df$index
stopifnot(all(uimg %in% all.ind))

jhut2.img = img


############## 
#####Trying to aggregate harvard's cortical and subcortical maps
##### into one unified atlas
######### not completed

# hoxall.img = hoxsubcort.img
# remove = hoxsubcort.df[
# 	grep("White Matter|Cortex", hoxsubcort.df$Label), 
# 	"index"]
# hoxall.img[ hoxall.img %in% remove] = 0

# subcortical = hoxall.img > 0               
# cortical = hoxcort.img > 0
# overlap = cortical & subcortical
# either= cortical | subcortical


# hoximg = hoxsubcort.img *0
# hoximg[(!overlap) & subcortical] = hoxsubcort.img[(!overlap) & subcortical]
# hoximg[(!overlap) & cortical] = hoxcort.img[(!overlap) & cortical]

# ####if it were ventricle, but overlapped, use cortical
# vents = hoxsubcort.df$index[grepl("Ventricle", hoxsubcort.df$Label)]
# hoximg[overlap & (hoxsubcort.img %in% vents)] = 
# 	hoxcort.img[overlap & (hoxsubcort.img %in% vents)]


# stopifnot(sum(overlap) == 0)
# # indices = hoxall.img[overlap]

# rm(list=c("subcortical", "cortical", "overlap"))




# hoxall.img = hoxall.img + (100*hoxcort.img)

# hoxall.df = cbind(hoxsubcort.df, study="Subcortical")
# hoxall.df = rbind(hoxall.df, cbind(hoxcort.df, study="Cortical"))

get.ind = function(img, df){
	ind.list = llply(df$index, function(x){
		return(which(img == x))
	}, .progress = "text")
	return(ind.list)
}

tal.list = get.ind(tal.img, tal.df)
mni.list = get.ind(mni.img, mni.df)
hoxcort.list = get.ind(hoxcort.img, hoxcort.df)
hoxsubcort.list = get.ind(hoxsubcort.img, hoxsubcort.df)

jhut1.list = get.ind(jhut1.img, jhut1.df)
jhut2.list = get.ind(jhut2.img, jhut2.df)


save(tal.df, tal.img, tal.list,
	mni.df, mni.img, mni.list,
	hoxcort.df, hoxcort.img, hoxcort.list,
	hoxsubcort.df, hoxsubcort.img, hoxsubcort.list,
	jhut2.df, jhut1.df, jhut2.img, jhut1.img,
	jhut1.list, jhut2.list,
	file= file.path(outdir, "All_FSL_Atlas_Labels.Rda"))



