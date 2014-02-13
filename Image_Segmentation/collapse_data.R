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

mat$ss.img <- paste0(mat$img, "_SS_First_Pass_Mask_0.1")
mat$ss.img <- file.path(dirname(mat$ss.img), 
  "Skull_Stripped", 
  basename(mat$ss.img))

set.seed(20131113)
## all combinations
nfiles <- nrow(mat)
combos <- combn(nrow(mat), floor(nrow(mat)/2))
subset <- sample(1:ncol(combos), 1)
ind <- 1:nfiles
train.ind <- ind %in% combos[, subset]
test.ind <- !train.ind

irow <- 1
rda <- paste0(mat$img[irow], "_Data.rda")
load(file=rda)
df$fname <- mat$img[irow]
# all.df <- df
all.df <- vector(mode="list", length=nfiles)
all.df[[1]] <- df
rm(list="df")
gc()
for (irow in 2:nfiles){
  rda <- paste0(mat$img[irow], "_Data.rda")
  load(file=rda) 
  df$fname <- mat$img[irow]
  # all.df <- rbind(all.df, df)
  all.df[[irow]] <- df
  rm(list="df")
  for (i in 1:10) gc()
  print(irow) 
}

all.df <- rbindlist(all.df)

df <- all.df

rda <- file.path(basedir, "Collaped_Data.rda")
save(list = c("df"), file=rda)

