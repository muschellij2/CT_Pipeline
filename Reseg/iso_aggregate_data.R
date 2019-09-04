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

##############################
# Keeping files where predictors exist
##############################

##############################
# Run lmod number of models - not all the models - leave out
##############################
fdf.run = fdf[ fdf$group == "Train", ]

L = nrow(fdf.run)
all.df = vector(
    mode = "list", 
    length = L)    
l.keep.ind = vector(mode="list", 
    length = L)
imod = 1

for (imod in seq(L)){
    img.pred = readRDS(fdf.run$outfile[imod])

    df = img.pred$df
    keep.ind = img.pred$keep.ind
    nim = img.pred$nim

    df = df[ keep.ind, ]
    med.ztemp = median(df$zscore_template)

    df$gr_medztemp = (df$zscore_template > 
        med.ztemp) 

    l.keep.ind[[imod]] = keep.ind

    df$img = fdf.run$img[imod]
    all.df[[imod]] = df
    rm(list=c("img.pred", "df"))
    print(imod)
}

all.df = bind_rows(all.df)

keep.colnames = colnames(all.df)
keep.colnames = keep.colnames[ 
    !keep.colnames %in% "gr_medztemp"
    ]

stopifnot(all(all.df$Y %in% c(0, 1)))
names(l.keep.ind) = fdf.run$img

runnames = colnames(all.df)
nosmooth = c("any_zero_neighbor",
        "thresh", "pct_zero_neighbor")
runnames = runnames[ !(runnames %in% 
    c("mask", "Y", "img", nosmooth))]


quants = plyr::dlply(all.df[, 
    c(runnames, "Y")], plyr::.(Y), 
    function(x) {
    r=t(colQuantiles(as.matrix(x[, runnames]), 
        probs = c(0, 0.001, 0.005, 0.01, 
            0.99, 0.995, 0.999, 1)))
    colnames(r) = runnames
    r
}, .progress = "text")

est.cutoffs = quants$`1`[ 
    c("0.5%", "99.5%"), ]

fname = file.path(resdir, 
   "iso_aggregate_data_cutoffs.rda")

multiplier_col = "zval2"

save(est.cutoffs, quants, 
    multiplier_col,
    file = fname)

keepnames = colnames(est.cutoffs)
include = rep(TRUE, length=nrow(all.df))
for (icut in keepnames){
    qcuts = est.cutoffs[, icut]
    colname = paste0(icut, ".cutoff")
    all.df[, colname] = 
        all.df[, icut] >= qcuts[1] & 
        all.df[, icut] <= qcuts[2]
    include = include & all.df[, colname]
    print(icut)
}



sum(all.df$Y[!include])/ sum(all.df$Y)
sum(include)/ length(include)
all.df$include.all = include

all.df$include = all.df$value >= 30 & 
    all.df$value <= 100


all.df$zval = all.df[, "zscore3.cutoff"] & 
    all.df$include &
    all.df$pct_thresh.cutoff
all.df$zval2 = all.df[, "zscore2.cutoff"] & 
    all.df$zval
all.df$zval_all = 
    all.df[, "zscore_template.cutoff"] & 
    all.df$zval2
all.df$zval2_medztemp = all.df$zval2 & 
    all.df$gr_medztemp 


sum(all.df$Y[! all.df$zval2]) / sum(all.df$Y)
sum(1-all.df$Y[! all.df$zval2]) / sum(1-all.df$Y)

all.df$multiplier = all.df[, 
    multiplier_col]


all.df$candidate = all.df$multiplier | all.df$Y == 1

seed = 20141022
set.seed(seed)
ich = which(all.df$Y == 1 & 
    all.df$multiplier) 
noich = which(all.df$Y != 1 & 
    all.df$multiplier)

size = 1e5
prop = .25
n.ich = ceiling(size*prop)
n.noich = size - n.ich
ich.ind = sample(ich, size=n.ich)
noich.ind = sample(noich, size=n.noich)
samp.ind = sort(c(ich.ind, noich.ind))

# samp.ind = sample(nrow(df), size= 1e4)
samps = seq(nrow(all.df)) %in% samp.ind

mult.df = all.df[, !colnames(all.df) %in% keep.colnames]
all.df = all.df[, colnames(all.df) %in% keep.colnames]    


img = all.df$img
all.df$img = NULL

# train = all.df[samps,]
# test = all.df[!samps,]

fname = file.path(resdir, 
    "iso_aggregate_data.rda")

save(all.df, mult.df, samps, img, 
    size, prop, n.ich, n.noich, 
    fdf.run, seed, l.keep.ind,
    file = fname)

train.img = img[samps]
train = all.df[samps,]
# tr.img = img[samps,]
train.mult.df = mult.df[samps,]

test.img = img[!samps]
test = all.df[!samps,]
test.mult.df = mult.df[!samps, ]

test = test[ test.mult.df$candidate, ]
test.img = test.img[test.mult.df$candidate]

test.mult.df = test.mult.df[ test.mult.df$candidate, ]

fname = file.path(resdir, 
    "iso_candidate_aggregate_data.rda")

save(train.mult.df, train, 
    test, test.mult.df,
    train.img, test.img,
    multiplier_col = multiplier_col,
    file = fname)    

