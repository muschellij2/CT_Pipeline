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
rm(list="train")

rda <- file.path(basedir, "Models", paste0("Aggregate_Model.rda"))
load(file=rda)

xmod.mat <- mod.mat
mods <- mod.mat
cn <- colnames(mods)
cn <- cn[ !(cn %in% "converged")]
mods <- mods[,cn]

nscen <- nrow(scenarios)
auc <- rep(NA, length=nscen)

### number of folds is number of cores
nfolds <- 75
folds <-c(sapply(1:nfolds, rep, length=ceiling(nscen/nfolds)))
folds <- folds[1:nscen]

# mods <- matrix(rnorm( (length(cols)+1) * nscen), 
# 	ncol=length(cols)+1, 
# 	nrow=nscen)
# colnames(mods) <- c("(Intercept)", cols)

cn <- colnames(mods)
cn <- cn[ !(cn %in% "converged")]
v <- as.data.frame(valid)
v$"(Intercept)" <- 1
v <- v[, cn]
v <- as.matrix(v)

### 
mods <- t(mods)
yy <- y <- valid$y

rm(list=c("valid", "scenarios"))



ifold <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if (is.na(ifold)) ifold <- 2

# for (ifold in unique(folds)){
	mod.ind <- folds == ifold
	mod.mat <- mods[, mod.ind]
	conv <- xmod.mat[mod.ind, "converged"]
	# mod.mat <- mod.mat[, 1:10]
	system.time(preds <- v %*% mod.mat)

	rm(list=c("v", "mod.mat", "mods"))


	# y <- list(y)
	fpr.stop <- .1
	auc <- matrix(nrow=ncol(preds), ncol=2)

	system.time({
		# i <- 1	
		# auc <- apply(preds, 2, function(pred) {
		for (icol in 1:ncol(preds)){
			# print(i)
			# i <<- i + 1		
			if (conv[icol] == 1){
				pred <- preds[, icol]
				rpred <- prediction(pred, y)
				rperf <- performance(rpred, "auc", fpr.stop=fpr.stop)
				pauc <- as.numeric(rperf@y.values)/fpr.stop
				rperf <- performance(rpred, "auc")
				full.auc <- as.numeric(rperf@y.values)
				auc[icol,] <- c(full.auc, pauc)
			} 
			print(icol)
		}
		# })
	})


# }


rda <- file.path(basedir, "Models", paste0("AUC_Fold_", ifold, ".rda"))
save(list=c("auc", "folds", "ifold"), file=rda)


