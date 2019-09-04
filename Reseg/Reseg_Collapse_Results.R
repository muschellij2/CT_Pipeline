####################################
## This code is for aggregate models
## CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
####################################
#####################################
rm(list=ls())
library(plyr)
library(dplyr)
library(cttools)
library(fslr)
library(ggplot2)
set.seed(20150518)
homedir = "/Applications"
rootdir = file.path("/Volumes/DATA_LOCAL", 
    "Image_Processing")
if (Sys.info()[["user"]] %in% 
    "jmuschel") {
  homedir = "~"
  rootdir = file.path(
    "/legacy/dexter/disk2/smart", 
    "stroke_ct", "ident")
}
progdir = file.path(rootdir, 
    "programs")
basedir = file.path(rootdir, 
    "Registration")
tempdir = file.path(rootdir, 
    "Template")
atlasdir = file.path(tempdir, 
    "atlases")
outdir = file.path(basedir, 
    "results")

segdir = file.path(progdir, 
    "Reseg")
source(file.path(segdir, 
    "Reseg_performance_functions.R"))
source(file.path(segdir, 
    "ggplot_smooth_func.R"))

correct = "Rigid"

options = c("Rigid")
# options = c("none", "Rigid")


#### load voxel data
outfile = file.path(outdir, 
    "Reseg_111_Filenames.Rda")
load(file = outfile)

gsubber = function(cname){
  fdf[, cname] = gsub("^/dexter", "/legacy/dexter", fdf[, cname])
  fdf
}
for (icol in colnames(fdf)) {
  fdf = gsubber(icol)
}

# for (correct in options){
if ("all.df" %in% ls()){
    rm(list="all.df")
}
for (i in 1:3) gc()
correct = match.arg(correct, options)
adder = switch(correct, 
    "none"= "",
    "Rigid" = "_Rigid")
fdf$roi.fname = switch(correct,
    "none" = fdf$roi,
    "Rigid" = fdf$rig_ssroi
    )
fdf$img.fname = switch(correct,
        "none" = fdf$img,
        "Rigid" = fdf$rig_ssimg   
        )
fdf$res_rda_stub = file.path(fdf$outdir, 
    paste0("Reseg_Results_", 
        nii.stub(fdf$img.fname, bn = TRUE)
        )
    )

fdf$vol_rda_stub = file.path(fdf$outdir, 
            paste0("Reseg_", 
                nii.stub(fdf$img, bn=TRUE), 
                "_estvolume")
            )

fdf$truevol_rda = file.path(fdf$outdir, 
        paste0("Reseg_", 
            nii.stub(fdf$img, bn=TRUE), 
            "_truevolume.Rda"))


filename = file.path(outdir, 
    paste0("Reseg_Result_Formats", 
        adder, ".Rda"))
load(filename)

makedir = sapply( fdf$outdir, 
    function(x) {
    if (!file.exists(x)){
        dir.create(x, showWarnings =FALSE)
    }
})

##############################
# Keeping files where predictors exist
##############################
fdf$outfile = nii.stub(basename(fdf$img))
fdf$outfile = paste0("Reseg_", fdf$outfile, 
    "_predictors", adder, ".Rda")
fdf$outfile = file.path(fdf$outdir, 
    fdf$outfile)
fdf = fdf[file.exists(fdf$outfile), ]

###########################################
# Get cutoffs and make multiplier
###########################################
mods = c("logistic", "lasso", 
    "rf", 
    "gam")
    # , 
    # "cforest")
# , "gam", "cforest")
imod = mods[1]



iimg <- as.numeric(
    Sys.getenv("SGE_TASK_ID"))
if (is.na(iimg)) iimg = 21

revver = function(x){
    ss = strsplit(x, "[.]")
    ss = lapply(ss, rev)
    ss = sapply(ss, paste, collapse = ".")
}

iapp = c("_native")

all.results = NULL
x = fdf[iimg,]

for (iimg in seq(nrow(fdf))) {
    x = fdf[iimg,]

    for (iapp in c("", "_native")){
        
        for (imod in mods){
            # stub_fname = file.path(outdir, 
            #     paste0(
            #      "Reseg_Aggregate_models",
            #      adder))
            # check_stub =  paste0(stub_fname,
            #     "_", imod)
            # check_fname = paste0(check_stub,
            #     ".Rda")   
            
            res_rda = paste0(x$res_rda_stub, 
                "_", imod, iapp, ".Rda")
            if (
                file.exists(res_rda) 
                ) {

                ########################
                # Use dcmtables to calc
                # Volume for native data
                ########################
                if (iapp == "_native"){  
                    vol_rda_stub = paste0(
                        x$vol_rda_stub, 
                        "_", imod, ".Rda")
                    load(vol_rda_stub)
                    est_vols = 
                        est_vols * 1000
                    truevol_rda = 
                        x$truevol_rda
                    xxx = load(truevol_rda)
                    true_vols = est_vols
                    true_vols[
                        seq_along(true_vols)
                        ] = 
                        truevol * 1000
                    names(est_vols) = paste0(
                        "estvol.", 
                        names(est_vols)) 
                    names(true_vols) = paste0(
                        "truevol.", 
                        names(true_vols)) 
                    all_vols = c(est_vols, 
                        true_vols)
                    rm(list=xxx) 
                }
                load(file = res_rda)
                xres = results
                results = data.frame(
                    t(unlist(results)))
                colnames(results) = revver(
                    names(results)
                    )
                results$mod = imod
                results$iimg = iimg
                results$app = iapp
                if (iapp == "_native") {  
                    cn = names(true_vols)
                    results[1, cn] = 
                        true_vols
                }
                all.results = rbind(
                    all.results, 
                    results)
            }   
            print(imod)
        }
        print(iapp)
    }


    print(iimg)
}


all.results = filter(all.results, 
    !mod %in% "cforest")


long = reshape(all.results, 
    direction = "long",
    idvar = c("iimg", "mod", "app"),
    varying = revver(names(unlist(xres))),
    timevar = "cutoff"
    )
long$id = fdf$id[long$iimg]
long$group = fdf$group[long$iimg]
rownames(long) = NULL

long$mod = factor(long$mod)
long$cutoff = factor(long$cutoff)
long$app[ long$app %in% "" ] = "Rigid"
long$app[ 
    long$app %in% "_native" ] = "Native"

long$tvol = long$truevol / 1000
long$evol = long$estvol / 1000

outrda = file.path(outdir, 
    paste0("Reseg_Results.Rda"))
save(long, 
    all.results, 
    fdf, 
    file = outrda)


