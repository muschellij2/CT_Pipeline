####################################
## This code is for aggregate models
## CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
####################################
#####################################
rm(list=ls())
library(fslr)
library(ichseg)
library(ROCR)
library(matrixStats)
library(randomForest)
library(methods)
library(parallel)
library(doMC)
set.seed(20150518)
cores <- 8

cl <- makeCluster(cores)
registerDoMC(cores=cores)
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

irow = 1
x = fdf[irow,]

##############################
# Keeping files where predictors exist
##############################
fname = file.path(resdir, 
    "iso_candidate_aggregate_data.rda")
x = load(fname)


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
mod = strip_model(mod)
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

stub_fname = file.path(resdir, 
    "iso_aggregate_models")
fname = paste0(stub_fname, 
    "_logistic.rda")

save(modlist, 
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
test.rf.pred = test.rf.pred * 
test$multiplier

rf.modlist = mod_func(
    test.rf.pred, 
    test$Y, 
    fpr.stop)
rf.modlist$mod = rf.mod
rf.modlist$mod.time = rf.time
rf.modlist$mod.perf = NULL

modlist = rf.modlist

fname = paste0(stub_fname, 
    "_rf.rda")

save(modlist,  
    file = fname)

