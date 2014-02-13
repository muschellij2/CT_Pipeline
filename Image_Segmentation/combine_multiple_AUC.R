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
source(file.path(progdir, "pred.prob.R"))


rda <- file.path(basedir, "First_Train_Data.rda")
load(file=rda)

rda <- file.path(basedir, "Models", paste0("AUC_Fold_1.rda"))
load(file=rda)

all.auc <- rep(NA, length=length(folds))

for (ifold in 1:max(folds)){
	rda <- file.path(basedir, "Models", paste0("AUC_Fold_", ifold, ".rda"))
	load(file=rda)
	mod.ind <- folds == ifold
	all.auc[mod.ind] <- auc
	print(ifold)
}
rm(list=c("folds", "ifold", "auc"))

top.ind <- which(all.auc > 0.98)
top.auc <- all.auc[top.ind]
top.ind <- top.ind[order(-top.auc)]
top.scen <- scenarios[top.ind,]
rowSums(top.scen)
colSums(top.scen)
	
	# print(iscen)
# # }
