#######################################
## This code is for aggregate prediction 
## of Image Segmentation of 
## CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
########################################
#####################################
## This code is for aggregate data 
## CT
## Author: John Muschelli
#################################
rm(list=ls())

library(methods)
library(neurobase)
library(dplyr)
library(matrixStats)
library(ichseg)
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  rootdir = Sys.getenv("dex")
}
basedir = file.path(rootdir, 
    "PITCH_reconverted", 
    "processed_data")
resdir = file.path(rootdir, 
    "PITCH_reconverted", 
    "results")
filename = file.path(resdir, 
    "iso_filename.rds")
fdf = readRDS(filename)

irow = 1
x = fdf[irow,]

groups = c("Validation", "Test")
igroup = as.numeric(
    Sys.getenv("SGE_TASK_ID"))
if (is.na(igroup)) {
    igroup = 1 
}
group = groups[igroup]
##############################
# Run lmod number of models - 
# not all the models - leave out
##############################
ffdf = fdf
fdf.run = ffdf[ffdf$group == group,]

L = nrow(fdf.run)
all.df = vector(
    mode = "list", 
    length = L)   

l.keep.ind = vector(mode="list", 
    length = L)
imod = 1

fname = file.path(resdir, 
   "iso_aggregate_data_cutoffs.rda")
load(fname)

for (imod in seq(L)){
    img.pred = readRDS(fdf.run$outfile[imod])

    df = img.pred$df
    keep.ind = img.pred$keep.ind
    nim = img.pred$nim

    df = df[ keep.ind, ]
    med.ztemp = median(df$zscore_template)

    df$gr_medztemp = (df$zscore_template > 
        med.ztemp) 

    keep.colnames = colnames(df)

    stopifnot(all(df$Y %in% c(0, 1)))

    keepnames = colnames(est.cutoffs)
    include = rep(TRUE, length=nrow(df))
    for (icut in keepnames){
        qcuts = est.cutoffs[, icut]
        colname = paste0(icut, ".cutoff")
        df[, colname] = 
            df[, icut] >= qcuts[1] & 
            df[, icut] <= qcuts[2]
        include = include & df[, colname]
    }

    sum(df$Y[!include])/ sum(df$Y)
    sum(include)/ length(include)
    df$include.all = include

    df$include = 
        df$value >= 30 & df$value <= 100


    df$zval = df[, "zscore3.cutoff"] & 
        df$include &
        df$pct_thresh.cutoff
    df$zval2 = df[, "zscore2.cutoff"] & 
        df$zval
    df$zval_all = 
        df[, "zscore_template.cutoff"] & 
        df$zval2
    df$zval2_medztemp = df$zval2 & 
        df$gr_medztemp 


    sum(df$Y[! df$zval2]) / sum(df$Y)
    sum(1-df$Y[! df$zval2]) / sum(1-df$Y)

    sum(df$Y[! df$zval2_medztemp]) / 
        sum(df$Y)
    sum(1-df$Y[! df$zval2_medztemp]) / 
        sum(1-df$Y)

    df$multiplier = df[, multiplier_col]
    # zval2

    df$candidate = 
        df$multiplier | df$Y == 1

    df = df[ df$candidate, ]

    set.seed(20150215)
    size = min(nrow(df), 1e5)
    samp = sample(nrow(df), size=size)
    df = df[samp, ]

    # l.keep.ind[[imod]] = keep.ind

    df$img = fdf.run$img[imod]
    df$mask = df$mask > 0
    all.df[[imod]] = df
    rm(list=c("img.pred", "df"))
    print(imod)
}

all.df = bind_rows(all.df)

# samp.ind = sample(nrow(df), size= 1e4)
mult.df = all.df[, 
    !colnames(all.df) %in% keep.colnames]
all.df = all.df[, 
    colnames(all.df) %in% keep.colnames]


img = all.df$img
all.df$img = NULL

fname = file.path(resdir, 
    paste0("iso_", tolower(group), 
        "_data.rda"))

save(all.df, mult.df, img, 
    fdf.run, 
    file = fname)

