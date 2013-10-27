################################
# Written 2013Oct24
# Author: John Muschelli
# Purpose: Create per-scan regression models with 10% of the data
# Output: Model rda, ROC Curve, Probability Map
# Use of output: intensity-based normalization of CT
################################
rm(list=ls())
library(oro.dicom)
library(plyr)
library(AnalyzeFMRI)
library(car)
library(ROCR)
basedir <- "/Volumes/Seagate Backup Plus Drive/Image_Processing/Test_5"
if (Sys.info()[["user"]] %in% "jmuschel") basedir <- "/dexter/disk2/smart/stroke_ct/ident/Test_5"
progdir <- file.path(dirname(basedir), "programs")
source(file.path(progdir, "Zscore.R"))

files <- list.files(path=basedir, pattern="*.nii.gz", full.names=TRUE)
stubs <- gsub(".nii.gz", "", files, fixed=TRUE)
mat <- matrix(stubs, ncol=2, byrow=TRUE)
mat <- data.frame(mat, stringsAsFactors=FALSE)

### make sure the data only has images and ROIS
colnames(mat) <- c("img", "roi")
test <- gsub("_ROI", "", mat$roi)
stopifnot(all(test == mat$img))


for (irow in 1:nrow(mat)){
  truth <- readNIfTI(mat$roi[irow], reorient=FALSE)
  
  img1 <- readNIfTI(mat$img[irow], reorient=FALSE)
  # system.time(zs <- zscore2(img1))
  system.time(zax <- zscore(img1))
  system.time(zcor <- zscore(img1, margin=2))
  system.time(zsag <- zscore(img1, margin=1))
  mask <- img1 >= 0 & img1 <= 100
  
  df <- data.frame(val=img1[mask],
                   zax=zax[mask], 
                   zcor=zcor[mask],
                   zsag=zsag[mask],
                   y = truth[mask])
  df$over40 <- df$val >= 40
  
  # over0 <- img1 > 0
  # system.time(smooth <- GaussSmoothArray(img1, mask=over0))
  system.time(smooth <- GaussSmoothArray(img1))
  df$smooth <- smooth[mask]
  
  perc <- 0.1
  N <- nrow(df)
  samp <- sample(1:N, floor(N*perc), replace=FALSE)
  
  train <- df[samp,]
  test <- df[ !(1:N %in% samp), ]
  mod <- glm(y ~ ., data=train, family=binomial)
  summary(mod)
  
  pred <- predict(mod, newdata=test, type="response")
  sum(train$y)
  
  rpred <- prediction(predictions=pred, labels=test$y)
  rperf <- performance(rpred, "tpr", "fpr")
  
  ### make ROC Curve
  outdir <- file.path(dirname(mat$img[irow]), "prob_maps")  
  fname <- basename(mat$img[irow])
  pdf(file.path(outdir, paste0(fname, ".pdf")))
    plot(rperf, main=fname)
  dev.off()
  
  ### save results to rda 
  mod$data <- NULL
  pimg <- img1
  pimg@.Data[!mask] <- 0
  all.pred <- predict(mod, newdata=df, type="response")
  pimg@.Data[mask] <- all.pred
  rda <- paste0(mat$img[irow], ".rda")
  save(list = c("mod", 
            "df", 
            "samp", 
            "perc"),
            file=rda)
  
  ### make a probability map (datatype=16 means float)
  pfname <- file.path(outdir, fname)
  rimg <- range(pimg)
  pimg@cal_max <- rimg[2]
  pimg@cal_min <- rimg[1]
  pimg@scl_inter <- 0
#   pimg@data_type <- "double"
  pimg2 <- nifti(pimg, datatype=16)
  pimg2@vox_offset
  writeNIfTI(pimg2, filename=pfname)
  print(irow)
  
}