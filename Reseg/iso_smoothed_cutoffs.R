####################################
## This code is for aggregate models
## CT
## Author: John Muschelli
####################################
####################################
rm(list=ls())
library(methods)
library(ichseg)
library(neurobase)
library(dplyr)
library(matrixStats)
library(randomForest)
library(readr)
library(ROCR)
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
    "iso_aggregate_data.rda")
load(fname)

cut_fname = file.path(resdir, 
   "iso_aggregate_data_cutoffs.rda")
x = load(cut_fname)

fdf = fdf %>% 
    mutate(
    pimg = file.path(outdir,
        paste0(mod_type, "_prob.nii.gz")),
    sm_pimg = sub("_prob", "_sm_prob", 
        pimg))

cutoff = modlist$mod.dice.coef[1,
    'cutoff']

fdf.run = fdf %>% 
    filter(group == "Train")

L = nrow(fdf.run)


all.p = vector(
    mode = "list",
    length = L
    ) 
iimg = 1
################################
# Read smoothed probability image
# concatenate
################################        
for (iimg in seq(L)){
    idf = fdf.run[iimg,]
    out_fname = idf$outfile
    img.pred = read_rds(out_fname)

    df = img.pred$df
    keep.ind = img.pred$keep.ind
    nim = img.pred$nim

    df$candidate = ich_candidate_voxels(df, 
        cutoffs = est.cutoffs)

    img = readnii(idf$sm_pimg)
    df$p = c(img)

    df = df[keep.ind,]

    all.p[[iimg]] = df %>% 
        select(p, candidate, Y)
    print(iimg)
}

all.p = bind_rows(all.p)

stopifnot(nrow(all.p) == nrow(all.df))

################################
# Take the subsampled data
# So comparable to unsmoothed
################################
test = all.p[ !samps, ]
test = test[ test$candidate, ]

################################
# Performance Metrics
################################
modlist = mod_func(test$p, 
    test$Y, 
    fpr.stop = 0.01)
modlist$mod.perf = NULL

stub_fname = file.path(resdir, 
    "iso_aggregate_models")
outname = paste0(stub_fname, 
    "_", mod_type, "_smoothed.rda")

save(modlist,  
    file = outname)

