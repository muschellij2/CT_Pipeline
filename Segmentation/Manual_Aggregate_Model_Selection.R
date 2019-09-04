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

fpr.stop = 0.01
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
            "thresh")
    runnames = runnames[ !(runnames %in% 
        c("mask", "Y", "img", nosmooth))]


    ########################
    # Standardizing Data
    ########################    
    cmeans = colMeans(train[, runnames])
    csds = apply(train[, runnames], 2, sd)


    train[, runnames] = t(t(train[, runnames])/csds - cmeans/csds)
    test[, runnames] = t(t(test[, runnames])/csds - cmeans/csds)

    all.mods = list()
    runmod = function(formstr, family= binomial()){
        form = as.formula(formstr)
        mod = glm(formula=form, data=train, family=family)
        return(mod)
    }
    form = formstr = "Y ~ . - mask"

    orig.mod = runmod(formstr) 
    orig.smod = summary(orig.mod)

    orig.qmod = runmod(formstr, family=quasibinomial()) 
    orig.sqmod = summary(orig.qmod)    


    all.mods = c(all.mods, orig.mod=list(orig.mod))

    pvals = coef(orig.smod)[, "Pr(>|z|)"]
    pvals = pvals[ names(pvals) != "(Intercept)"]

    v = vif(orig.mod)
    keep = names(v)[v <= 5]    

    formstr = paste0("Y ~", paste0(keep, collapse = " + "))

    vifmod = runmod(formstr) 
    svifmod = summary(vifmod)

    all.mods = c(all.mods, vifmod=list(vifmod))

    orig.mod = remove_lmparts(orig.mod)
    vifmod = remove_lmparts(vifmod)

    #####################################
    # Run all subsets regression
    #####################################
    # vars = names(pvals)

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
    # # all.mods = c(all.mods, list(mod))


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
    # all.mods = c(all.mods, list(mod))



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
    # all.mods = c(all.mods, list(mod))


    # formstr = paste0("Y ~ dist_centroid", 
    #     "+ smooth10 + smooth20 + zscore_template")
    # #########################################
    # # Run model without high VIF
    # #########################################
    # mod = runmod(formstr)
    # smod = summary(mod)
    # smod
    # vif(mod)

    # all.mods = c(all.mods, list(mod))

    # formstr = paste0("Y ~ dist_centroid", 
    #     "+ smooth10 + zscore_template")
    # #########################################
    # # Run model without high VIF
    # #########################################
    # mod = runmod(formstr)
    # smod = summary(mod)
    # smod
    # vif(mod)
    # all.mods = c(all.mods, list(mod))

    reses = array(dim=c(3,3, length(all.mods)))
    for (imod in seq_along(all.mods)){
        runmod = all.mods[[imod]]

        ############################
        # train predictions
        ############################
        train.pred = predict(runmod, train, type="response")
        train.pred = train.pred * train.mult.df$multiplier        
        tr.pred = prediction(train.pred, train$Y)

        dice.cutoff = get_max_dice(obj=tr.pred)[, "cutoff"]

        acc.cutoff = get_acc(tr.pred)[, 'cutoff']

        sens.cutoff = get_senscut(obj = tr.pred, 
            fpr.stop = fpr.stop)[, "cutoff"]

        ############################
        # train predictions
        ############################
        test.pred = predict(runmod, test, type="response")
        test.pred = test.pred * test.mult.df$multiplier

        # t.pred = prediction(test.pred, test$Y)
        cut = c(dice.cutoff, acc.cutoff, sens.cutoff)
        tabs = lapply(cut, function(icut){
                extrantsr:::my.tab(test.pred > icut, test$Y)
        })

        res = matrix(ncol=length(tabs), nrow=3)
        rownames (res) = colnames(res)= c("dice", "acc", "sens")
        res[, "dice"] = sapply(tabs, get.dice)
        res[, "acc"] = sapply(tabs, get.dice)
        res[, "sens"] = sapply(tabs, get.sens)

        reses[,,imod ] = res
    
    }
    dimnames(reses) = list(cutoff=c("dice", "acc", "sens"), 
        maesure=c("dice", "acc", "sens"), model=NULL)

    apply(reses, 3, function(res){
        apply(res, 2, which.max)
    })    

    maxes = apply(reses, 3, colMaxs)
    rownames(maxes) = c("dice", "acc", "sens")
    # test$orig.pred = test$orig.pred * test.mult.df$multiplier
    # opred = prediction(test$orig.pred, test$Y)

    # train$orig.pred = predict(orig.mod, train, type="response")
    # train$orig.pred = train$orig.pred * train.mult.df$multiplier

    # tr.opred = prediction(train$orig.pred, train$Y)

    # #########################################
    # # Run model without high VIF
    # #########################################
    # train$pred = predict(mod, train, type="response")
    # train$pred = train$pred * train.mult.df$multiplier

    # tr.pred = prediction(train$pred, train$Y)

    # test$pred = predict(mod, test, type="response")
    # test$pred = test$pred * test.mult.df$multiplier

    # pred = prediction(test$pred, test$Y)


    # ################# 
    # # Training data
    # #################
    # get_max_dice(obj=tr.pred)
    # get_max_dice(obj=tr.opred)

    # get_acc(tr.pred)
    # get_acc(tr.opred)

    # get_pauc(tr.pred, fpr.stop = fpr.stop)
    # get_pauc(tr.opred, fpr.stop = fpr.stop)

    # get_senscut(obj = tr.pred, fpr.stop = fpr.stop)

    # get_senscut(obj = tr.opred, fpr.stop = fpr.stop)

    # cuts = lapply(list(tr.pred, tr.opred), function(x){
    #     dice.cutoff = get_max_dice(obj=x)[, "cutoff"]

    #     acc.cutoff = get_acc(x)[, 'cutoff']

    #     sens.cutoff = get_senscut(obj = x, 
    #         fpr.stop = fpr.stop)[, "cutoff"]

    #     c(dice.cutoff = dice.cutoff, 
    #         acc.cutoff = acc.cutoff,
    #         sens.cutoff = sens.cutoff)
    # })


    # ###########################
    # # Run cutoffs on test data
    # ###########################
    # tabs = mapply(function(tpred, cut){
    #     lapply(cut, function(icut){
    #         extrantsr:::my.tab(tpred > icut, test$Y)
    #     })
    # }, list(test$pred, test$orig.pred), cuts)

    # dices = matrix(ncol=ncol(tabs), nrow=nrow(tabs))
    # senss = accs = dices
    # for (icol in seq(ncol(tabs))){
    #     dices[, icol] = sapply(tabs[,icol], get.dice)
    #     accs[, icol] = sapply(tabs[,icol], get.acc)
    #     senss[, icol] = sapply(tabs[,icol], get.sens)
    # }




    ################# 
    # Test data
    #################



    # get_max_dice(obj=pred)
    # get_max_dice(obj=opred)

    # get_acc(pred)
    # get_acc(opred)

    # get_pauc(pred, fpr.stop = fpr.stop)
    # get_pauc(opred, fpr.stop = fpr.stop)

    # get_senscut(obj = pred, fpr.stop = fpr.stop)

    # get_senscut(obj = opred, fpr.stop = fpr.stop)
