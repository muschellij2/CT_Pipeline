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
library(cttools)
library(fslr)
library(ROCR)
library(matrixStats)
library(mgcv)
library(extrantsr)
library(randomForest)
library(methods)
library(glmnet)
library(partykit)
library(parallel)
library(doMC)
set.seed(20150518)
cores <- 8
cl <- makeCluster(cores)
registerDoMC(cores=cores)
homedir = "/Applications"
rootdir = file.path("/Volumes/DATA_LOCAL", 
    "Image_Processing")
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = file.path("/legacy/dexter/disk2/smart", 
    "stroke_ct", "ident")
}
progdir = file.path(rootdir, "programs")
basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")
outdir = file.path(basedir, "results")

segdir = file.path(progdir, "Reseg")
source(file.path(segdir, 
    "Reseg_performance_functions.R"))

correct = "Rigid"

options = c("Rigid")
# options = c("none", "Rigid")


#### load voxel data
outfile = file.path(outdir, 
    "Reseg_111_Filenames.Rda")
load(file = outfile)

# for (correct in options){
if ("all.df" %in% ls()){
    rm(list="all.df")
}
for (i in 1:3) gc()
correct = match.arg(correct, options)
adder = switch(correct, 
    "none"= "",
    "Rigid" = "_Rigid")


filename = file.path(outdir, 
    paste0("Reseg_Result_Formats", 
        adder, ".Rda"))
load(filename)

makedir = sapply( fdf$outdir, function(x) {
    if (!file.exists(x)){
        dir.create(x, showWarnings =FALSE)
    }
})
irow = 1
x = fdf[irow,]

##############################
# Keeping files where predictors exist
##############################
outfiles = nii.stub(basename(fdf$img))
outfiles = paste0("Reseg_", outfiles, 
    "_predictors", adder, ".Rda")
outfiles = file.path(fdf$outdir, 
    outfiles)
fdf = fdf[file.exists(outfiles), ]


fname = file.path(outdir, 
    paste0(
    "Reseg_Candidate_Aggregate_data", 
        adder, ".Rda"))
load(fname)


train$gr_medztemp = NULL
train$mode = fdf$mode[match(train.img, 
    fdf$img)]
test$mode = fdf$mode[match(test.img, 
    fdf$img)]

runnames = colnames(train)
nosmooth = c("any_zero_neighbor",
        "thresh", "pct_zero_neighbor")
runnames = runnames[ !(runnames %in% 
    c("mask", "Y", "img", nosmooth))]


runmod = function(formstr){
    form = as.formula(formstr)
    mod = glm(formula=form, 
        data=train, 
        family=binomial())
    return(mod)
}
sds = colSds(as.matrix(train))
names(sds) = colnames(sds)
novar = sds == 0
novar = colnames(train[novar])
novar = novar[!novar %in% c("mask")]
novar = c("", novar)
novar = paste0(novar, collapse=" - ")

keepers = c("smooth10", "zscore3", 
    "pct_thresh", "value", "dist_centroid", 
    "smooth5", "flipped_value", "prob_img", 
    "win_z", "moment1", "zscore_template")
keepers = paste0(keepers, collapse = " + ")
formstr = paste0("Y ~ ", keepers)

test$multiplier = test.mult.df[, 
    "multiplier"]


mod.time = system.time({
    mod = runmod(formstr)
})
smod = summary(mod)
mm = model.matrix(mod)
mm = mm[ ,
    !colnames(mm) %in% c("(Intercept)")
    ]

N = nrow(test)
fpr.stop = 0.01

######################
# Standard logistic Model
######################
mod = keep_mod(mod)
test.pred = predict(mod, 
    test, 
    type = "response")
test.pred = test.pred * test$multiplier

###################################
# Get measures and cutoffs
###################################
modlist = mod_func(
    test.pred, 
    test$Y, 
    fpr.stop)
modlist$mod = mod
modlist$mod.time = mod.time
modlist$mod.perf = NULL

mod.modlist = modlist

stub_fname = file.path(outdir, 
    paste0("Reseg_Aggregate_models", 
        adder))
fname = paste0(stub_fname, 
    "_logistic_reduced.Rda")

save(fdf.run, 
    modlist, 
    smod, 
    file = fname)

##############################
# Random Forest Model
##############################
rf.time <- system.time({
    rf.mod = randomForest(x = mm,
    y = factor(train$Y), 
    do.trace = TRUE)
    })

test.rf.pred = rep(0, 
    length=nrow(test))
test.rf.pred[test$multiplier] = 
    as.numeric(
    predict(rf.mod, test[test$multiplier,], 
        type="prob")[,"1"]
)

##############################
# Random Forest Predictions
##############################
cat("# randomForest Prediction \n")

test.rf.pred = as.numeric(test.rf.pred)
test.rf.pred[ test.rf.pred > 1] = 1
test.rf.pred = test.rf.pred * test$multiplier

rf.modlist = mod_func(
    test.rf.pred, 
    test$Y, 
    fpr.stop)
rf.modlist$mod = rf.mod
rf.modlist$mod.time = rf.time
rf.modlist$mod.perf = NULL

modlist = rf.modlist

fname = paste0(stub_fname, 
    "_rf_reduced.Rda")

save(fdf.run, 
    modlist,  
    file = fname)

