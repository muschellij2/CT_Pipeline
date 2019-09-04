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
library(foreach)
library(itertools)
library(doParallel)
set.seed(20150518)
cores <- 8
cl <- makeCluster(cores)
registerDoParallel(cl)
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

subtest = test[test$multiplier,]

num_splits = 80
predictions <- foreach(d = 
    isplitRows(subtest, chunks=num_splits),
          .combine=c, 
          .packages=c("partykit")) %dopar% {
    tt = predict(cforest.mod, newdata=d)
    tt = do.call("rbind", tt)[, 2]
    tt
}


test.cforest.pred[ test$multiplier ] = 
    predictions
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

stub_fname = file.path(outdir, 
    paste0("Reseg_Aggregate_models", 
        adder))
fname = paste0(stub_fname, 
    "_cforest.Rda")

save(fdf.run, 
    modlist,  
    file = fname)

# all.df$img = img

# }
if (!is.null(cl)) stopCluster(cl)
