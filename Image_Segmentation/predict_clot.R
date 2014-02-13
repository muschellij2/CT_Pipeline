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
library(ggplot2)
basedir <- "/Volumes/DATA_LOCAL/Image_Processing/Test_5"
if (Sys.info()[["user"]] %in% "jmuschel") basedir <- "/dexter/disk2/smart/stroke_ct/ident/Test_5"
rootdir <- dirname(basedir)
progdir <- file.path(rootdir, "programs")
source(file.path(progdir, "Zscore.R"))
source(file.path(progdir, "List_Moment_Generator.R"))
source(file.path(progdir, "pred.prob.R"))
source(file.path(progdir, "fslhd.R"))

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
### test to see if all files exist
imgs <- unlist(mat)
imgs <- paste0(imgs, ".nii.gz")
stopifnot(all(file.exists(imgs)))

# irow <- grep("102323", mat$img)
irow <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if (is.na(irow)) irow <- 1

### test to see if all files exist
check <- mat[irow,]
check <- paste0(check, ".nii.gz")
stopifnot(all(file.exists(check)))

# for (irow in 1:nrow(mat)){
  ## manual segmentation
  outdir <- file.path(dirname(mat$img[irow]), "prob_maps")  

  img1 <- readNIfTI(mat$img[irow], reorient=FALSE)
  fname <- basename(mat$img[irow])

  ss.mask <- readNIfTI(mat$ss.img[irow], reorient=FALSE)
  mask <- ss.mask > 0

  rda <- paste0(mat$img[irow], "_Data.rda")
  load(rda)

  set.seed(20131113)
  perc <- 0.1
  N <- nrow(df)
  samp <- sample(1:N, floor(N*perc), replace=FALSE)
  
  train <- df[samp,]
  test <- df[ !(1:N %in% samp), ]
  print(sum(train$y))
  mod <- glm(y ~ val + zax + 
    izax + izcor + izsag + 
    mean + sd + skew +
    over40 +
    smooth1 + smooth3 + smooth5 + smooth10 + smooth20,
    data=train, family=binomial)
  summary(mod)

  step.mod <- stepAIC(mod)

  # mod2 <- glm(y ~ . - zcor -zsag, 
  #   data=train, family=binomial)

  # mod3 <- glm(y ~ ., 
  #   data=train, family=binomial)

  # mod4 <- glm(y ~ . - zcor - zsag, 
  #   data=train, family=binomial)

  # mod5 <- glm(y ~ . - mean - sd - skew, 
  #   data=train, family=binomial)

  # mod2 <- glm(y ~ . - val, data=train, family=binomial)
  # summary(mod2)
  fpr <- 0.1
  get.pred <- function(model, auc=FALSE, fpr.stop = 0.1){
    pred <- predict(model, newdata=test, type="response")
    # sum(train$y)
  
    rpred <- prediction(predictions=pred, labels=test$y)
    if (auc) {
      rperf <- performance(rpred, "auc", fpr.stop=fpr.stop)
    } else {
      rperf <- performance(rpred, "tpr", "fpr", fpr.stop=fpr.stop)
    }
    return(rperf)
  }
  # mods <- list(mod, mod2, mod3, mod4, mod5, step.mod)
  mods <- list(mod, step.mod)
  mods <- lapply(mods, function(mod) {
    scrape.mod(mod)
  })
  # preds <- lapply(mods, function(model) 
    # predict(model, newdata=test, type="response") )
  preds <- lapply(mods, function(model) {
    pred.prob(model, test) 
  })
  rpreds <- lapply(preds, function(pred) 
    prediction(predictions=pred, labels=test$y) )  
  rperfs <-  lapply(rpreds, function(rpred) 
    performance(rpred, "tpr", "fpr", fpr.stop=fpr) ) 

  pAUC <- function(rpred, fpr) {
    auc <- performance(rpred, "auc", fpr.stop=fpr)
    auc <- auc@y.values[[1]]
    auc <- auc / fpr
  }
  fprs <- c(0.05, 0.1)
  aucs <-  sapply(fprs, function(fpr) 
      sapply(rpreds, pAUC, fpr=fpr) )
  colnames(aucs) <- fprs

  ### make ROC Curve
  pdf(file.path(outdir, paste0(fname, ".pdf")))
  plot(rperfs[[1]])
  for (imod in seq_along(rperfs)){
    plot(rperfs[[imod]], add=TRUE, col=imod)
  }
  for (imod in seq_along(rperfs)){
    plot(rperfs[[imod]], main=fname)
    abline(v=0.05, col="red")
    abline(v=0.1, col="blue")
    abline(h=0.9, col="black")
    text(x=0.5, y=0.5, 
      paste0("AUC at 0.05: ", 
        round(aucs[1, "0.05"], 3)), 
      col="red")
    text(x=0.5, y=0.4, 
      paste0("AUC at 0.1: ", 
        round(aucs[2, "0.1"], 3)), 
      col="blue") 
    vars <- attr(mods[[imod]]$terms, "term.labels")  
    text(.8, seq(.8, 0.2, length=length(vars)), vars)
  }
  dev.off()
  
  rda <- paste0(mat$img[irow], "_Models.rda")
  save(list = c("mods", 
            "samp", 
            "perc"),
            file=rda)

  ### save results to rda  - clear data to free up space
  # mod$data <- NULL
  all.pred <- pred.prob(mod, test=df)
  pimg <- img1
  pimg@.Data[!mask] <- 0
  #all.pred <- array(all.pred, dim=dim(pimg@.Data))
  pimg@.Data[mask] <- all.pred
  # pimg@.Data <- all.pred

  ### make a probability map (datatype=16 means float)
  pfname <- file.path(outdir, fname)
  rimg <- range(pimg)
  pimg@cal_max <- rimg[2]
  pimg@cal_min <- rimg[1]
  pimg@scl_inter <- 0
  pimg@datatype <- 16
  pimg@bitpix <- 32
#   pimg@data_type <- "double"
  # pimg2 <- nifti(pimg, datatype=16)
  # pimg2@vox_offset
  writeNIfTI(pimg, filename=pfname)
  print(irow)
  
# }