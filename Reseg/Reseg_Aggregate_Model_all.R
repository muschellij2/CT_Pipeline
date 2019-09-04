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
  rootdir = file.path("/dexter/disk2/smart", 
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
formstr = paste0("Y ~ . - mask", novar) 

test$multiplier = test.mult.df[, 
    "multiplier"]


mod.time = system.time({
    mod = runmod(formstr)
})
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
    "_logistic.Rda")

save(fdf.run, 
    modlist,  
    file = fname)

###################################
# LASSO Model
###################################
lasso.time = system.time({
    lasso.mod <- cv.glmnet(mm, train$Y, 
    type.measure = "class", 
    family = "binomial",
    parallel = TRUE)
    })

test.mm = model.matrix(
    as.formula(formstr), 
    test)
test.mm = test.mm[, colnames(mm)]

test.lasso.pred = predict(lasso.mod, 
    newx = test.mm,
    s = "lambda.1se", 
    type = "response")[, "1"]
test.lasso.pred[!test$multiplier] = 0
names(test.lasso.pred) = NULL
test.lasso.pred[ test.lasso.pred > 1] = 1

lasso.modlist = mod_func(
    test.lasso.pred, 
    test$Y, 
    fpr.stop)
lasso.modlist$mod = lasso.mod
lasso.modlist$mod.time = lasso.time
lasso.modlist$mod.perf = NULL

modlist = lasso.modlist

fname = paste0(stub_fname, 
    "_lasso.Rda")

save(fdf.run, 
    modlist,  
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
    "_rf.Rda")

save(fdf.run, 
    modlist,  
    file = fname)


### used gam - but bam supposedly faster
gam.time = system.time({
    gam.mod = bam(Y ~ 
    s(moment1) + 
    s(moment2) + 
    s(skew) + 
    s(kurtosis) + 
    s(value) + 
    thresh +
    s(zscore1) + 
    s(zscore2) + 
    s(zscore3) + 
    s(pct_thresh) + 
    pct_zero_neighbor + 
    any_zero_neighbor +
    s(dist_centroid) +
    s(smooth10) +
    s(smooth20) 
    , data=train, family= binomial(), 
    method = "fREML")
})
# + mode
gam.time


##############################
# GAM PREDs
##############################
test.gam.pred = rep(0, length=nrow(test))
test.gam.pred[test$multiplier] = as.numeric(
    predict(gam.mod, 
        test[test$multiplier,], 
        type="response")
)

cat("GAM Prediction \n")

test.gam.pred = as.numeric(test.gam.pred)
test.gam.pred[ test.gam.pred > 1] = 1
test.gam.pred = test.gam.pred * test$multiplier

gam.modlist = mod_func(
    test.gam.pred, 
    test$Y, 
    fpr.stop)
gam.modlist$mod = gam.mod
gam.modlist$mod.time = gam.time
gam.modlist$mod.perf = NULL

modlist = gam.modlist

fname = paste0(stub_fname, 
    "_gam.Rda")

save(fdf.run, 
    modlist,  
    file = fname)


###########################
### Conditional Forest
###########################
cforest.time <- system.time({
     cforest.mod = cforest(
        factor(Y) ~ . - mask, 
       data = train,
       cores = cores)
     })    
cforest.time 
###################################
# cforest Model
###################################
test.cforest.pred = rep(0, 
    length=nrow(test))

tt = 
    predict(cforest.mod, 
        newdata = test[test$multiplier,],
        type = "prob")
tt = do.call("rbind", tt)[, 2]
test.cforest.pred[ test$multiplier ] = tt
test.cforest.pred[!test$multiplier] = 0
names(test.cforest.pred) = NULL
test.cforest.pred[ 
    test.cforest.pred > 1 ] = 1

cforest.modlist = mod_func(
    test.cforest.pred, 
    test$Y, 
    fpr.stop)
cforest.modlist$mod = cforest.mod
cforest.modlist$mod.time = cforest.time
cforest.modlist$mod.perf = NULL    

modlist = cforest.modlist

fname = paste0(stub_fname, 
    "_cforest.Rda")

save(fdf.run, 
    modlist,  
    file = fname)

# all.df$img = img

fname = paste0(stub_fname, ".Rda")

save(
    # mods, 
    fdf.run, 

    mod.modlist,
    gam.modlist,
    cforest.modlist,
    lasso.modlist,
    rf.modlist
    file = fname)
# }
