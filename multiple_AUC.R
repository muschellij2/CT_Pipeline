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

rda <- file.path(basedir, "Models", paste0("Aggregate_Model.rda"))
load(file=rda)

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

v <- as.data.frame(valid)
v$"(Intercept)" <- 1
v <- v[, colnames(mods)]
v <- as.matrix(v)

### 
mods <- t(mods)
y <- valid$y

y1 <- y == 1
y0 <- !y1

ifold <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if (is.na(ifold)) ifold <- 1

# for (ifold in unique(folds)){
	mod.ind <- folds == ifold
	mod.mat <- mods[, mod.ind]
	preds <- v %*% mod.mat


	# system.time({
	# 	i <- 1	
	# 	auc <- apply(preds, 2, function(pred) {
	# 		print(i)
	# 		i <<- i + 1
	# 		colAUC(pred, y, alg="ROC")
	# 	})	
	# })

	fpr.stop <- .1
	system.time({
		i <- 1	
		auc <- apply(preds, 2, function(pred) {
			print(i)
			i <<- i + 1			
			rpred <- prediction(pred, y)
			rperf <- performance(rpred, "auc", fpr.stop=fpr.stop)
			pauc <- as.numeric(rperf@y.values)/fpr.stop
			rperf <- performance(rpred, "auc")
			full.auc <- as.numeric(rperf@y.values)
			c(full.auc, pauc)
		})	
	})


# }


rda <- file.path(basedir, "Models", paste0("AUC_Fold_", ifold, ".rda"))
save(list=c("auc", "folds", "ifold"), file=rda)


	
	# print(iscen)
# # }
