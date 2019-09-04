###################################################################
## This code is for Image Segmentation of CT
## The code is R but calls out FSL
##
## Author: John Muschelli
## Last updated: May 20, 2014
###################################################################
###################################################################
rm(list=ls())
library(plyr)
library(cttools)
library(devtools)
library(ROCR)
library(fslr)
library(mgcv)
library(extrantsr)
library(getopt)
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

segdir = file.path(progdir, "Segmentation")
source(file.path(segdir, "performance_functions.R"))


my.tab <- function(
  x, 
  y, 
  dnames=c("x", "y")) {
  x = as.numeric(x)
  y = as.numeric(y)
  stopifnot(all(unique(c(x,y)) %in% c(0, 1, NA)))
  tt = sum(x * y)
  t1=sum(x)
  t2=sum(y)
  tab = matrix(c(length(x)-t1-t2+tt,  t1-tt, t2-tt, tt), 2, 2)
  n = list(c("FALSE", "TRUE"), c("FALSE", "TRUE"))
  names(n) = dnames
  dimnames(tab) = n
  tab = as.table(tab)
  return(tab) 
}

correct = "Rigid"
# options = c("none", "N3", "N4", "N3_SS", "N4_SS",
#         "SyN", "SyN_sinc", "Rigid", "Affine", "Rigid_sinc", 
#         "Affine_sinc")
options = c("none", "N3_SS", "N4_SS", 
		"Rigid", "Rigid_sinc")

spec = matrix(c(
	'correct', 'c', 1, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

correct = opt$correct
print(opt)

# options = c("Rigid", "Rigid_sinc")

#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )

outfile = file.path(outdir, "111_Filenames_with_volumes.Rda")
load(file = outfile)


# iscen <- as.numeric(Sys.getenv("SGE_TASK_ID"))
# if (is.na(iscen)) iscen = 1


# scenarios = expand.grid(iimg = seq(nrow(fdf)), 
# 	correct = options, stringsAsFactors = FALSE )
# iimg = scenarios$iimg[iscen]
# correct = scenarios$correct[iscen]

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


iimg <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(iimg)) iimg = 1

# for (correct in options){

runx = x = fdf[iimg,]

fname= nii.stub(basename(x$img))
fname = paste0(fname, "_predictors", adder, "_training.Rda")
outfile = file.path(x$outdir, fname)
dsets = load(outfile)

N = nrow(fdf)
n.ind = 1e4
all.train = data.frame(
	matrix(NA, nrow=n.ind * N, ncol = ncol(train)+3)
	)
colnames(all.train) = c(colnames(train), "id", "pred", "spred")

rm(list=dsets)

filename = file.path(outdir, 
    paste0("Result_Formats", adder, ".Rda"))
load(filename)

mod.outdir = fdf.run$outdir[lmod]
moddname = nii.stub(basename(fdf.run$img[lmod]))
moddname = file.path(mod.outdir, 
	paste0(moddname, "_models", adder, ".Rda"))
xx = load(moddname)
mod = mods$mod 

ind = 0
fpr.stop = 0.01
for (iimg in seq(N)){

	runx = x = fdf[iimg,]
	    

	iddir = x$iddir
	id.outdir = x$outdir	    
	print(iimg)

	###################################
	# Load Predictors
	###################################
	fname= nii.stub(basename(x$img))
	fname = paste0(fname, "_predictors", adder, "_training.Rda")
	outfile = file.path(x$outdir, fname)

	xxx = load(file=outfile)

	train$id = x$id
	train$pred = short_predict(mod, train)


	##################
	# Getting Smoothed prediction
	##################
	cn = "mod_agg"
	type= "_zval2"
	outimg = nii.stub(x$img, bn=TRUE)
	outimg = file.path(id.outdir, 
		paste0(outimg, "_", cn, adder))	

	outimg = paste0(outimg, "_smoothed", adder, type)
	spred = readNIfTI(outimg, reorient=FALSE )

	train$spred = spred@.Data[keep.ind[samps]]

	end.ind = ind + n.ind
	all.train[ seq(ind + 1, end.ind), ] = train

	ind = end.ind

	rm(list=xxx)
	# train$truevol = x$truevol
	# train$varslice = x$varslice
}


drop.ids = fdf.run$id
drop.ids = drop.ids[!is.na(drop.ids)]


nopred = which(fdf$id %in% drop.ids)
ffdf = fdf[-nopred, ]
nr = nrow(ffdf)
valid.ind = ceiling(nr/2)
test.ind = seq( valid.ind +1, nr)
valid.ind = seq(1, valid.ind)

group = "Training"
if (group == "Training"){
	subset.ind = valid.ind
}
if (group == "Test"){
	subset.ind = test.ind
}
varslice = ffdf$varslice[subset.ind]
subset.ids = ffdf$id[subset.ind]
all.truevol =  ffdf$truevol
truevol = ffdf$truevol[subset.ind]


cut.filename = file.path(outdir, 
paste0("Model_Cutoffs", adder, ".Rda"))

xcuts = load(file=cut.filename)

cuts = sapply(xcuts, function(obj){
	obj = get(obj)
	obj['mod_agg']
} )

scut.filename = file.path(outdir, 
paste0("Smooth_Model_Cutoffs", adder, ".Rda"))

ycuts = load(file=scut.filename)

scuts = sapply(ycuts, function(obj){
	obj = get(obj)
	obj['mod_agg']
} )

##################################
# Subset the data to sets
##################################
train = all.train[ all.train$id %in% subset.ids, ]

pred = prediction(train$pred, train$Y)
perf_dice = dice(pred)

##################################
# Use the pROC curve
##################################
perf <- limit_pauc(pred, fpr.stop)

##################################
# Get results for prediction
##################################
N = nrow(train)
sens.cut = get_senscut(pred, fpr.stop=fpr.stop, 
    N = N,
    predictions=train$pred, Y = train$Y)
    print(sens.cut)
pauc = get_pauc(pred, fpr.stop=fpr.stop)
    print(pauc)
dice.coef = get_max_dice(pred)
    print(dice.coef)
acc = get_acc(pred)
    print(acc)


##################################
# Get smoothed results for prediction
##################################
spred = prediction(train$spred, train$Y)
sperf_dice = dice(spred)


ssens.cut = get_senscut(spred, fpr.stop=fpr.stop, 
    N = N,
    predictions=train$spred, Y = train$Y)
    print(ssens.cut)
spauc = get_pauc(spred, fpr.stop=fpr.stop)
    print(spauc)
sdice.coef = get_max_dice(spred)
    print(sdice.coef)
sacc = get_acc(spred)
    print(sacc)




##################################
# Run measures on the cutoffs
##################################
cut.tabs = lapply(cuts, function(x) {
	my.tab(train$pred > x, train$Y)
})

cut.stabs = lapply(scuts, function(x) {
	my.tab(train$spred > x, train$Y)
})
	