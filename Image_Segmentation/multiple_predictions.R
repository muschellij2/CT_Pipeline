################################
# Written 2013Oct24
# Author: John Muschelli
# Purpose: Extract models and make predictions
# Output: ROC Curve
# Use of output: prediction model for classifier
################################

rm(list=ls())
library(data.table)
# library(pROC)
basedir <- "/Volumes/DATA_LOCAL/Image_Processing/Test_5"
if (Sys.info()[["user"]] %in% "jmuschel") basedir <- "/dexter/disk2/smart/stroke_ct/ident/Test_5"
progdir <- file.path(dirname(basedir), "programs")
# source(file.path(progdir, "Zscore.R"))
# source(file.path(progdir, "List_Moment_Generator.R"))
source(file.path(progdir, "pred.prob.R"))

# # test.ind <- !train.ind

# train <- df[samp,]

rda <- file.path(basedir, "First_Train_Data.rda")
load(file=rda)


# test <- df[ !(1:nrow(df) %in% samp), ]

# ifold <- as.numeric(Sys.getenv("SGE_TASK_ID"))

# if (is.na(ifold)) ifold <- 1

nscen <- nrow(scenarios)
### number of folds is number of cores
nfolds <- 75
folds <-c(sapply(1:nfolds, rep, length=ceiling(nscen/nfolds)))
folds <- folds[1:nscen]

train <- as.data.frame(train)


ifold <- 1
mod.mat <- matrix(0, ncol=length(cols)+2, nrow=nscen)
colnames(mod.mat) <- c("(Intercept)", cols, "converged")

################################################################
####### Loop over folds - extract coefficients #################
######## and put into a matrix and then save   #################
################################################################

for (ifold in 1:nfolds){
	rda <- file.path(basedir, "Models", 
		paste0("Model_Fold_", ifold, ".rda"))
	load(file=rda)
	# coef(mod)

	scen.ind <- which(folds == ifold)
	for (iscen in scen.ind){
		coefs <- mods[[iscen]]$coefficients
		mod.mat[iscen, names(coefs)] <- coefs
		mod.mat[iscen, "converged"] <- mods[[iscen]]$converged
		# print(iscen)
	}
	print(ifold)
}


rda <- file.path(basedir, "Models", paste0("Aggregate_Model.rda"))
save(list="mod.mat", file=rda)


