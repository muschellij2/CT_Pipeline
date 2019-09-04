###################################################################
## This code is for aggregate prediction of Image Segmentation of 
## CT, with all subsets
##
## Author: John Muschelli
## Last updated: May 20, 2014
###################################################################
###################################################################
rm(list=ls())
library(plyr)
library(cttools)
library(fslr)
library(ROCR)
library(matrixStats)
library(mgcv)
library(extrantsr)
library(car)
library(getopt)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")
outdir = file.path(basedir, "results")
rerun = FALSE

segdir = file.path(progdir, "Segmentation")
source(file.path(segdir, "performance_functions.R"))

spec = matrix(c(
    'correct', 'c', 1, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

correct = opt$correct
print(opt)

options = c("none", "N3_SS", "N4_SS", 
      "Rigid", "Rigid_sinc")

if (is.null(correct)) correct = "Rigid"

# for (correct in options){
    if ("all.df" %in% ls()){
        rm(list="all.df")
    }
    for (i in 1:3) gc()
    correct = match.arg(correct, options)
    adder = switch(correct, 
        "none"= "",
        "N3"="_N3",
        "N4" = "_N4",
        "N3_SS" = "_N3_SS",
        "N4_SS" = "_N4_SS", 
        "SyN" = "_SyN",
        "SyN_sinc" = "_SyN_sinc",
        "Rigid" = "_Rigid",
        "Affine" = "_Affine",
        "Rigid_sinc" = "_Rigid_sinc",
        "Affine_sinc" = "_Affine_sinc")


    ichunk <- as.numeric(Sys.getenv("SGE_TASK_ID"))
    if (is.na(ichunk)) ichunk = 1

    outfile = file.path(outdir, 
            paste0("All_Subsets_Aggregate_Model_Results_Chunk_",
            ichunk, 
                adder, ".Rda"))
if (!file.exists(outfile) | rerun ){
    #################################
    # Load candidate data
    #################################
    fname = file.path(outdir, 
        paste0("Candidate_Aggregate_data", adder, ".Rda"))
    xx = load(fname)
    rm(list=c("train", "train.mult.df"))
    gc()

    fname = file.path(outdir, 
            paste0("All_Subsets_Aggregate_Model", 
                adder, ".Rda")) 
    load(fname)

    test$"(Intercept)" = 1
    Y = test$Y
    test = test[, colnames(coef.mat)]
    test = as.matrix(test)
    family = binomial()
    multiplier = test.mult.df$multiplier
    rm(list="test.mult.df")
    gc()



    ################################
    # Making chunks for looping
    ################################
    chunksize = 300
    lthresh=  .Machine$double.eps^0.5
    nchunks = ceiling(n.forms/chunksize)
    chunks = c(sapply(seq(1, nchunks), 
        rep, length=chunksize))
    chunks = chunks[seq(n.forms)]


    ################################
    # Making required summands
    ################################
    N = length(Y)
    n.pos = sum(Y)
    n.neg = N - n.pos
    CS = cumsum(rep(1, N))

    all.dice = all.acc = matrix(nrow=n.forms, ncol=2)
    all.pauc =  matrix(nrow=n.forms, ncol=1)
    all.senscut = matrix(nrow=n.forms, ncol=4)
    fpr.stop = 0.01


    # for (ichunk in seq(nchunks)){
        ind = which(chunks == ichunk)
        mat = as.matrix(coef.mat[ind, ])
        preds = tcrossprod(test, mat)

        ################################
        # Looping over results
        ################################
        n = ncol(preds)
        ipred = 1
        x = preds[,ipred]
        x = family$linkinv(x)
        x[x < lthresh ] = 0
        xx = get_meas(x)
        dice = matrix(nrow = n, ncol = NCOL(xx$dice))
        acc = matrix(nrow = n, ncol = NCOL(xx$acc))
        pauc = matrix(nrow = n, ncol = NCOL(xx$pauc))
        senscut = matrix(nrow = n, ncol = NCOL(xx$senscut))
        pb = txtProgressBar(max=n, style=3)
        for (ipred in seq(n)){
            x = preds[,ipred]
            x = family$linkinv(x)
            x[x < lthresh ] = 0
            xx = get_meas(x)
            dice[ipred, ] = xx$dice
            acc[ipred, ] = xx$acc
            pauc[ipred, ] = xx$pauc
            senscut[ipred, ] = xx$senscut
            setTxtProgressBar(pb, ipred)
            # print(icol)
        }
        close(pb)       
        
        # pred.results <- aaply(preds, 2, function(x){
        #     x = family$linkinv(x)
        #     x[x < lthresh ] = 0
        #     get_meas(x)
        #     }, 
        #     .progress = "text")


        # dice = do.call("rbind", pred.results[, "dice"])
        # acc = do.call("rbind", pred.results[, "acc"])
        # pauc = do.call("rbind", pred.results[, "pauc"])
        # senscut = do.call("rbind", pred.results[, "senscut"])


        save(dice, acc, pauc, senscut, ind,
            file = outfile)
}
 