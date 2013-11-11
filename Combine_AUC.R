################################
# Written 2013Oct24
# Author: John Muschelli
# Purpose: Extract models and make predictions
# Output: ROC Curve
# Use of output: prediction model for classifier
################################

rm(list=ls())
library(data.table)
library(ROCR)
library(caTools)
# library(pROC)
basedir <- "/Volumes/DATA_LOCAL/Image_Processing/Test_5"
if (Sys.info()[["user"]] %in% "jmuschel") basedir <- "/dexter/disk2/smart/stroke_ct/ident/Test_5"
progdir <- file.path(dirname(basedir), "programs")
# source(file.path(progdir, "pred.prob.R"))
# source(file.path(progdir, "myperf.R"))

rda <- file.path(basedir, "First_Train_Data.rda")
load(file=rda)
vals <- train$val[train$y == 1]
rm(list="train")

nscen <- nrow(scenarios)

### number of folds is number of cores
nfolds <- 75
folds <-c(sapply(1:nfolds, rep, length=ceiling(nscen/nfolds)))
folds <- folds[1:nscen]
aucs <- NULL


for (imod in 1:nfolds){
	rda <- file.path(basedir, "Models", 
		paste0("AUC_Fold_", imod, ".rda"))
	load(file=rda)
	aucs <- rbind(aucs, auc)
	print(imod)
}

ranks <- rank(-aucs[,2])
models <- which(ranks <= 100)

top.mods <- scenarios[models,]
cn <- colnames(scenarios)
# cols <- apply(cn, )


best <- cbind(top.mods, pAUC=aucs[models,2], full.auc=aucs[models,1])
colSums(best[, cols])
