################################
# Written 2013Oct24
# Author: John Muschelli
# Purpose: Collapse data for the 5 cases and run a multi-subject model
# Output: ROC Curve
# Use of output: prediction model for classifier
################################

rm(list=ls())
library(oro.dicom)
library(plyr)
library(AnalyzeFMRI)
library(car)
library(ROCR)
library(ggplot2)
library(data.table)
library(biglm)
library(randomForest)
library(rpart)
basedir <- "/Volumes/DATA_LOCAL/Image_Processing/Test_5"
if (Sys.info()[["user"]] %in% "jmuschel") basedir <- "/dexter/disk2/smart/stroke_ct/ident/Test_5"
progdir <- file.path(dirname(basedir), "programs")
source(file.path(progdir, "Zscore.R"))
source(file.path(progdir, "List_Moment_Generator.R"))
source(file.path(progdir, "pred.prob.R"))


files <- list.files(path=basedir, pattern="*.nii.gz", full.names=TRUE)
stubs <- gsub(".nii.gz", "", files, fixed=TRUE)
stubs <- stubs[!grepl("Zero", stubs)]
stubs <- stubs[grepl("ROI", stubs)]
mat <- matrix(stubs, ncol=1, byrow=FALSE)
mat <- data.frame(mat, stringsAsFactors=FALSE)

mat[, 2] <- gsub("_ROI", "", mat[,1])
### make sure the data only has images and ROIS
colnames(mat) <- c("roi", "img")
mat <- mat[, c("img", "roi")]
test <- gsub("_ROI", "", mat$roi)
stopifnot(all(test == mat$img))
test <- NULL

mat$ss.img <- paste0(mat$img, "_SS_First_Pass_Mask_0.1")
mat$ss.img <- file.path(dirname(mat$ss.img), 
  "Skull_Stripped", 
  basename(mat$ss.img))


# ## all combinations
set.seed(20131106)
nfiles <- nrow(mat)
combos <- combn(nrow(mat), floor(nrow(mat)/2))
subset <- sample(1:ncol(combos), 1)
ind <- 1:nfiles
train.ind <- ind %in% combos[, subset]
test.ind <- !train.ind

all.mods <- list()
for (irow in 1:nrow(mat)){
  rda <- paste0(mat$img[irow], "_Models.rda")
  load(file=rda)
  all.mods <- c(all.mods, mods)
  print(irow)
}
formulas <- sapply(all.mods, formula)
uforms <- unique(formulas)

# run.mods <- mods

rda <- file.path(basedir, "Collaped_Data.rda")
load(file=rda)

  perc <- 0.01
  indices <- 1:nrow(df)
  dropscan <- mat$img[test.ind]

  drop.ind <- which(df$fname %in% dropscan )
  indices <- indices[ -drop.ind ]

  N <- length(indices)
  all.train <- df[indices,]
  samp <- sample(indices, floor(N*perc), replace=FALSE)
  
  train <- df[samp,]
  test <- df[ !(1:nrow(df) %in% samp), ]

  rm(list="df")

  train$y <- factor(train$y)

  rp <- rpart(formula=uforms[[1]], data=train)


  ntest <- nrow(test)
  nfolds <- 5
  folds <-c(sapply(1:nfolds, rep, length=ceiling(ntest/nfolds)))
  folds <- folds[1:ntest]

  # test$fold <- folds
  test$fname <- NULL

  test$pred <- NA
  ifold <- 1
  for (ifold in 1:nfolds){
    print(ifold)
    ind <- which(folds == ifold)
    tt <- test[ ind,  ]
    pred <- predict(rp, newdata=tt, type="class")
    test$pred[ind] <- pred
  }

  tab <- table(test$y, test$pred)
  print(prop.table(tab, 2))
  print(prop.table(tab, 1))

  rf <- randomForest(formula=uforms[[1]], data=train)
  rda <- file.path(basedir, "Rpart_and_RF.rda")
  save(rf, rp, file=rda)