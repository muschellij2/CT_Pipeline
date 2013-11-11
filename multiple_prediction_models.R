################################
# Written 2013Oct24
# Author: John Muschelli
# Purpose: Collapse data for the 5 cases and run a multi-subject model
# Output: ROC Curve
# Use of output: prediction model for classifier
################################

rm(list=ls())
library(data.table)
library(biglm)
basedir <- "/Volumes/DATA_LOCAL/Image_Processing/Test_5"
if (Sys.info()[["user"]] %in% "jmuschel") basedir <- "/dexter/disk2/smart/stroke_ct/ident/Test_5"
progdir <- file.path(dirname(basedir), "programs")
source(file.path(progdir, "Zscore.R"))
source(file.path(progdir, "List_Moment_Generator.R"))
source(file.path(progdir, "pred.prob.R"))

files <- list.files(path=basedir, pattern="*.nii.gz", full.names=TRUE)
stubs <- gsub(".nii.gz", "", files, fixed=TRUE)
mat <- matrix(stubs, ncol=2, byrow=TRUE)
mat <- data.frame(mat, stringsAsFactors=FALSE)

### make sure the data only has images and ROIS
colnames(mat) <- c("img", "roi")
test <- gsub("_ROI", "", mat$roi)
stopifnot(all(test == mat$img))

mat$ss.img <- paste0(mat$img, "_SS_First_Pass_Mask_0.1")
mat$ss.img <- file.path(dirname(mat$ss.img), 
  "Skull_Stripped", 
  basename(mat$ss.img))

### test to see if all files exist
imgs <- unlist(mat)
imgs <- paste0(imgs, ".nii.gz")
stopifnot(all(file.exists(imgs)))


rda <- file.path(basedir, "First_Train_Data.rda")
load(file=rda)

valid <- NULL

nscen <- nrow(scenarios)
auc <- rep(NA, length=nscen)

ifold <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if (is.na(ifold)) ifold <- 75

### number of folds is number of cores
nfolds <- 75
folds <-c(sapply(1:nfolds, rep, length=ceiling(nscen/nfolds)))
folds <- folds[1:nscen]

train <- as.data.frame(train)


nscen <- nrow(scenarios)

scen.ind <- which(folds == ifold)
# all.mod <- 	mod <- glm(y ~ ., 
# 		data=train[, c(cols, "y"), drop=FALSE], 
# 		family=binomial)
# coefs <- coef(all.mod)
mods <- vector(mode="list", length=nscen)
iscen <- scen.ind[860]

for (iscen in scen.ind){
	## get the variables for the model
	scen <- scenarios[iscen, ]
	vnames <- colnames(scen)[as.logical(scen)]
	vars <- cols[ cols %in% vnames ]
	# start.coef <- coefs[c("(Intercept)", vars)]
	vv <- vars
	vars <- c("1", vv)
	form <- paste0("y ~ ", paste0(vars, collapse="+"))
	form <- as.formula(form)
	system.time(mod <- bigglm(formula=form, 
		data=train, 
		family=binomial(),
		sandwich=FALSE, maxit=25))	
	# system.time(mod <- glm(formula=form, 
	# 	data=train, 
	# 	family=binomial()))
	mod <- scrape.mod(mod)
	# mod$vars <- c("(Intercept)", vv)
	mods[[iscen]] <- mod
	print(iscen)
}


rda <- file.path(basedir, "Models", paste0("Model_Fold_", ifold, ".rda"))
save(list = c("seed", "mods", "cols"), file=rda)

