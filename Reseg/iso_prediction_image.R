#####################################
## This code is for predicting
## Author: John Muschelli
#################################
rm(list=ls())
library(methods)
library(ichseg)
library(neurobase)
library(dplyr)
library(matrixStats)
library(randomForest)
library(extrantsr)
rootdir = file.path("/Volumes/DATA_LOCAL", 
    "Image_Processing")
if (Sys.info()[["user"]] %in% "jmuschel") {
  rootdir = Sys.getenv("dex")
}
basedir = file.path(rootdir, 
    "PITCH_reconverted", 
    "processed_data")
resdir = file.path(rootdir, 
    "PITCH_reconverted", 
    "results")
progdir = file.path(rootdir,
    "PITCH_reconverted",
    "code")
source(file.path(progdir, 
    "iso_performance_functions.R"))

filename = file.path(resdir, 
    "iso_filename.rds")
fdf = readRDS(filename)

mod_type = "rf"

stub_fname = file.path(resdir, 
    "iso_aggregate_models")
fname = paste0(stub_fname, 
    "_", mod_type, ".Rda")
load(fname)

cutoff = modlist$mod.dice.coef[1,
    'cutoff']
rm(modlist)
fname = paste0(stub_fname, 
    "_", mod_type, "_smoothed.rda")    
xx = load(fname)

scutoff = modlist$mod.dice.coef[1,
    'cutoff']

fdf = fdf %>% 
    mutate(
    pimg = file.path(outdir,
        paste0(mod_type, "_prob.nii.gz")),
    sm_pimg = sub("_prob", "_sm_prob", 
        pimg),
    pred_img = sub("prob[.]", "pred.",
        pimg),
    sm_pred_img = sub("prob[.]", "_sm_pred.",
        pimg),    
    )


##############################
# Keeping files where predictors exist
##############################
imod = as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(imod)) {
    imod = 1
}
L = nrow(fdf)

print(imod)
idf = fdf[imod,]

img = idf$img
mask_fname = sub("[.]nii",
    "_mask.nii", img)
if (!file.exists(mask_fname)) {
    ss = CT_Skull_Strip_robust(
        img, retimg = TRUE)
    mask = ss > 0
    writenii(mask, mask_fname)
} else {
    mask = readnii(mask_fname)
}
img = check_nifti(img)
ss = mask_img(img, mask)

template.file = system.file(
    "scct_unsmooth_SS_0.01.nii.gz", 
    package = "ichseg")

res = registration(
    filename = template.file, 
    template.file = ss, 
    retimg = TRUE, 
    typeofTransform = "Rigid", 
    reproducible = TRUE, 
    interpolator = "Linear", 
    remove.warp = FALSE
    )
get_native = function(in_img) {
    ants_apply_transforms(
    fixed = ss, moving = in_img, 
    interpolator = "genericLabel", 
    transformlist = res$fwdtransforms)
}

fnames = c(idf$pred_img, 
    idf$sm_pred_img)
fnames = c(fnames, sub("[.]nii",
    "_native.nii", fnames))
if (!all(file.exists(fnames))) {
    
    pimg = readnii(idf$pimg)     
    pred = pimg > cutoff
    cc = ants_bwlabel(img = pred, 
        k = 100, binary = TRUE)
    writenii(cc, idf$pred_img)
    cc = get_native(cc)
    writenii(cc, sub("[.]nii",
    "_native.nii", idf$pred_img))

    pimg = readnii(idf$sm_pimg)     
    pred = pimg > scutoff
    cc = ants_bwlabel(img = pred, 
        k = 100, binary = TRUE)
    writenii(cc, idf$sm_pred_img)
    cc = get_native(cc)
    writenii(cc, sub("[.]nii",
        "_native.nii", idf$sm_pred_img)) 
}

