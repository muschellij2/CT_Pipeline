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
library(corrplot)
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

if (is.null(correct)) correct = "none"

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

    #################################
    # Load candidate data
    #################################
    fname = file.path(outdir, 
        paste0("Candidate_Aggregate_data", adder, ".Rda"))
    xx = load(fname)
    
    runnames = colnames(train)
    nosmooth = c("any_zero_neighbor",
            "thresh", "pct_zero_neighbor")
    runnames = runnames[ !(runnames %in% 
        c("mask", "Y", "img", nosmooth))]

    runmod = function(formstr){
        form = as.formula(formstr)
        mod = glm(formula=form, data=train, family=binomial())
        return(mod)
    }
    form = formstr = "Y ~ . - mask"
    sds = colSds(as.matrix(train))
    names(sds) = colnames(sds)
    novar = sds == 0
    novar = colnames(train[novar])
    novar = novar[!novar %in% c("mask")]
    novar = c("", novar)
    novar = paste0(novar, collapse=" - ")
    formstr = paste0("Y ~ . - mask", novar) 

    orig.mod = runmod(formstr) 

    drop_int = function(x){
        x = x[x != "(Intercept)"]
    }

    cmat = cor(train[, drop_int(names(coef(orig.mod)))])

    corrplot(cmat, type = "upper")    
    cmat[lower.tri(cmat, diag=TRUE)] = NA
    smod = summary(orig.mod)

    # test$orig.pred = predict(orig.mod, test, type="response")
    # opred = prediction(test$orig.pred, test$Y)
    orig.mod = remove_lmparts(orig.mod)   

    pvals = coef(smod)[, "Pr(>|z|)"]
    pvals = pvals[ names(pvals) != "(Intercept)"]

    #####################################
    # Run all subsets regression
    #####################################
    vars = names(pvals)
    l = lapply(vars, function(x) c(FALSE, TRUE))
    names(l) = vars

    all.scen = expand.grid(l)
    forms = apply(all.scen, 1, function(x){
        s = sum(x)
        if (s == 0){
            xx = c("Y ~ 1")
        } else {
            v = vars[x]
            xx = paste0("Y ~ ", paste(v, collapse = "+ "))
        }
        return(xx)
    })
    n.forms = length(forms)
    coef.mat = all.scen
    coef.mat[!is.na(coef.mat)] = 0
    coef.mat = cbind("(Intercept)" = 0, coef.mat)
    # forms = forms[1:50]
    for (imod in seq_along(forms)){
        x = forms[imod]
        mod = runmod(x)
        cc = coef(mod)
        coef.mat[imod, names(cc)] = cc
        print(imod)
    }

    save(forms, coef.mat, n.forms,
        file = 
        file.path(outdir, 
            paste0("All_Subsets_Aggregate_Model", 
                adder, ".Rda")))    

    # test$"(Intercept)" = 1
    # Y = test$Y
    # test = test[, colnames(coef.mat)]
    # test = as.matrix(test)
    # family = binomial()




    # ################################
    # # Making chunks for looping
    # ################################
    # chunksize = 100
    # lthresh=  .Machine$double.eps^0.5
    # nchunks = ceiling(n.forms/chunksize)
    # chunks = c(sapply(seq(1, nchunks), 
    #     rep, length=chunksize))
    # chunks = chunks[seq(n.forms)]


    # ################################
    # # Making required summands
    # ################################
    # N = length(Y)
    # n.pos = sum(Y)
    # n.neg = N - n.pos
    # CS = cumsum(rep(1, N))
    # ichunk = 1

    # all.dice = all.acc = matrix(nrow=n.forms, ncol=2)
    # all.pauc =  matrix(nrow=n.forms, ncol=1)
    # all.senscut = matrix(nrow=n.forms, ncol=4)

    # for (ichunk in seq(nchunks)){
    #     ind = which(chunks == ichunk)
    #     mat = as.matrix(coef.mat[ind, ])
    #     preds = tcrossprod(test, mat)
    #     preds <- family$linkinv(preds)
    #     preds[ preds < lthresh] = 0


    #     pred.results <- aaply(preds, 2, get_meas, 
    #         .progress = "text")      


    #     dice = do.call("rbind", pred.results[, "dice"])
    #     acc = do.call("rbind", pred.results[, "acc"])
    #     pauc = do.call("rbind", pred.results[, "pauc"])
    #     senscut = do.call("rbind", pred.results[, "senscut"])

    #     all.dice[ind, ]= dice
    #     all.acc[ind, ]= acc
    #     all.pauc[ind, ]= pauc
    #     all.senscut[ind, ]= senscut
    #     print(ichunk)
    # }

    # colnames(all.dice) = colnames(dice)
    # colnames(all.acc) = colnames(acc)
    # colnames(all.senscut) = colnames(senscut)
    # colnames(all.pauc) = "pauc"



    # save(forms, coef.mat, 
    #     all.dice, all.acc, all.senscut, all.pauc,
    #     file = 
    #     file.path(outdir, 
    #         paste0("All_Subsets_Aggregate_Model", 
    #             adder, ".Rda")))

    # #########################################
    # # Taking out NAs
    # #########################################
    # formstr = paste0(form, 
    #     "- pct_zero_neighbor - any_zero_neighbor")

    # #########################################
    # # Run model without NAs
    # #########################################
    # mod = runmod(formstr)
    # smod = summary(mod)
    # vif(mod)

    # #########################################
    # # Get corelation matrix
    # #########################################
    # mm = model.matrix(as.formula(paste0(formstr, "-1")), 
    # data=train)
    # tcor = cor(mm)
    # diag(tcor) = NA
    # rm(list="mm")

    # high.ind = which(abs(tcor) > .85, arr.ind=TRUE)
    # high.cor = tcor[high.ind]
    # high.ind = data.frame(
    #     apply(high.ind, 2, function(x) colnames(tcor)[x])
    #     , stringsAsFactors=FALSE)
    # high.ind$cor = high.cor

    # formstr = paste0(form, 
    #     "- pct_zero_neighbor - any_zero_neighbor",
    #     "- moment4 - value - moment1 - smooth20")
    # #########################################
    # # Run model without high VIF
    # #########################################
    # mod = runmod(formstr)
    # smod = summary(mod)
    # smod
    # vif(mod)


    # formstr = paste0(form, 
    #     "- pct_zero_neighbor - any_zero_neighbor",
    #     "- moment4 - value - smooth20",
    #     "-zscore1", "-zscore2")
    # #########################################
    # # Run model without high VIF
    # #########################################
    # mod = runmod(formstr)
    # smod = summary(mod)
    # smod
    # vif(mod)


    # test$pred = predict(mod, test, type="response")

    # pred = prediction(test$pred, test$Y)


    # get_max_dice(obj=pred)
    # get_max_dice(obj=opred)

    # get_acc(pred)
    # get_acc(opred)

    # get_pauc(pred, fpr.stop = fpr.stop)
    # get_pauc(opred, fpr.stop = fpr.stop)

    # get_senscut(obj = pred, fpr.stop = fpr.stop)

    # get_senscut(obj = opred, fpr.stop = fpr.stop)
