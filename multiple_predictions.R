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

# irow <- grep("102323", mat$img)
iscen <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if (is.na(iscen)) iscen <- 1

nscen <- nrow(scenarios)
mods <- matrix(0, ncol=length(cols)+1, nrow=nscen)
colnames(mods) <- c("(Intercept)", cols)

for (iscen in 1:nscen){
# system.time({
	rda <- file.path(basedir, "Models", paste0("Model", iscen, ".rda"))
	load(file=rda)
	# coef(mod)

	var.classes <- get.stuff(mod)


	### checking to see if we used factors or not in the model
	stopifnot(all(var.classes %in% c('numeric', 'logical')))
	coefs <- coef(mod)
	have.log <- var.classes %in% "logical"
	if (any(have.log)){
	vars <- names(var.classes)[have.log]
	for (ivar in vars){
	  names(coefs) <- gsub(paste0(ivar, "TRUE"), ivar, names(coefs))
	}
	}

	mods[iscen, names(coefs)] <- coefs
	print(iscen)
}


	rda <- file.path(basedir, "Models", paste0("Aggregate_Model.rda"))
	save(list="mods", file=rda)
	# mods[[iscen]] <- mod
# })


# 	system.time(preds <- pred.prob(mod=mod, test=valid))
# 	# system.time(p2 <- predict(mod, newdata=valid))
# 	rpreds <- prediction(preds, valid$y)
# 	perf <- performance(rpreds, "auc")
# 	auc <- as.numeric(perf@y.values)
# 	print(iscen)
# # }
# 	rda <- file.path(basedir, "Models", paste0("AUC", iscen, ".rda"))
# 	save(list = c("auc", "preds"), file=rda)
# 	