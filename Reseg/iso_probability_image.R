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

mod = modlist$mod
fname = file.path(resdir, 
   "iso_aggregate_data_cutoffs.rda")
load(fname)

fdf = fdf %>% 
    mutate(
    pimg = file.path(outdir,
        paste0(mod_type, "_prob.nii.gz")),
    sm_pimg = sub("_prob", "_sm_prob", 
        pimg))

cutoff = modlist$mod.dice.coef[1,
    'cutoff']


##############################
# Keeping files where predictors exist
##############################
imod = as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(imod)) {
    imod = 1
}
L = nrow(fdf)

# for (imod in seq(L)) {
    print(imod)
    idf = fdf[imod,]

    fnames = c(idf$pimg, 
        idf$sm_pimg)
    if (!all(file.exists(fnames))) {
        img.pred = readRDS(idf$outfile)

        df = img.pred$df
        keep.ind = img.pred$keep.ind
        nim = img.pred$nim

        rm(img.pred)
        # df = df[ keep.ind, ]

        df$multiplier = 
        ich_candidate_voxels(
            df, 
            cutoffs = est.cutoffs)

        f = paste0("predict_", mod_type)
        func = get(f)
        p = func(mod, df[ df$multiplier, ])
        pimg = remake_img(p, nim, 
            df$multiplier)

        out_fname =  idf$pimg
        writenii(pimg, filename = out_fname)

        mask = remake_img(df$mask, nim)
        pimg = mask_img(pimg, mask)

        sm.pimg  = mean_image(pimg, 
            nvoxels = 1)
        sm.pimg[abs(sm.pimg) < 
            .Machine$double.eps^0.5 ] = 0
        sm.pimg = niftiarr(nim, sm.pimg)
        sm.pimg[is.na(sm.pimg)]= 0

        sout_fname = idf$sm_pimg
        writenii(sm.pimg, 
            filename = sout_fname)
        rm(df)
        rm(nim)
    }
# }
