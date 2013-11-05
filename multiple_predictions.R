################################
# Written 2013Oct24
# Author: John Muschelli
# Purpose: Extract models and make predictions
# Output: ROC Curve
# Use of output: prediction model for classifier
################################

rm(list=ls())
library(data.table)
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

mat$ss.img <- paste0(mat$img, "_SS_No1024_Mask_0.1")
mat$ss.img <- file.path(dirname(mat$ss.img), 
  "Skull_Stripped", 
  basename(mat$ss.img))

### test to see if all files exist
imgs <- unlist(mat)
imgs <- paste0(imgs, ".nii.gz")
stopifnot(all(file.exists(imgs)))


## all combinations
nfiles <- nrow(mat)
combos <- combn(nrow(mat), floor(nrow(mat)/2))
subset <- sample(1:ncol(combos), 1)
ind <- 1:nfiles

rda <- file.path(basedir, "Collaped_Data.rda")
load(file=rda)

rda <- file.path(basedir, "Collaped_Model.rda")
load(file=rda)

# test.ind <- !train.ind

# train <- df[samp,]
test <- df[ !(1:nrow(df) %in% samp), ]

seed <- 20131105
set.seed(seed)
## create valid/train data set
train.indices <- which(train.ind)
ntrain <- length(train.indices)
valid.ind <- sample(train.indices, floor(ntrain/2))

valid.fnames <- mat$img[valid.ind]
valid.ind <- which( train$fname %in% valid.fnames )

valid <- train[ valid.ind, ]
train <- train[ -valid.ind, ]

train <- as.data.frame(train)
valid <- as.data.frame(valid)

### get all combos of variables
cols <- colnames(train)[!colnames(train) %in% c("y", "fname")]
ncols <- length(cols)
scen <- lapply(1:ncols, function(x) c(FALSE, TRUE))
scenarios <- expand.grid(scen)
# irow <- grep("102323", mat$img)
iscen <- 1

nscen <- nrow(scenarios)
mods <- vector(length=nscen, mode="list")
for (iscen in 1:nscen){
	rda <- file.path(basedir, "Models", paste0("Model", iscen, ".rda"))
	load(file=rda)
	mods[[iscen]] <- mod
	preds <- pred.prob(mod=mod, test=test)
}
