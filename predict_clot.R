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
basedir <- "/Volumes/DATA_LOCAL/Image_Processing/Test_5"
if (Sys.info()[["user"]] %in% "jmuschel") basedir <- "/dexter/disk2/smart/stroke_ct/ident/Test_5"
progdir <- file.path(dirname(basedir), "programs")
source(file.path(progdir, "Zscore.R"))
source(file.path(progdir, "List_Moment_Generator.R"))

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

irow <- grep("101307", mat$img)

irow <- as.numeric(Sys.getenv("SGE_TASK_ID"))

# for (irow in 1:nrow(mat)){
  ## manual segmentation
  truth <- readNIfTI(mat$roi[irow], reorient=FALSE)
  
  img1 <- readNIfTI(mat$img[irow], reorient=FALSE)
  # system.time(zs <- zscore2(img1))
  system.time(izax <- zscore(img1))
  system.time(izcor <- zscore(img1, margin=2))
  system.time(izsag <- zscore(img1, margin=1))


  ss.mask <- readNIfTI(mat$ss.img[irow], reorient=FALSE)
  mask <- ss.mask > 0

  zero100 <- img1 >= 0 & img1 <= 100
### masking out non-brain 
  ss.img <- img1
  ss.img[!mask] <- 0

  system.time(mom <- image_moment(img1, mask=zero100, 
      nvox= 1, moment=c(1, 2, 3)))
  mom <- lapply(mom, function(x) {
    ind <- which(is.na(x), arr.ind=TRUE)
    x[ind] <- 0
    x
  })
### these averages are only with brain tissue - not overall
  system.time(zax <- zscore(ss.img))
  system.time(zcor <- zscore(ss.img, margin=2))
  system.time(zsag <- zscore(ss.img, margin=1))  
  


  df <- data.frame(val=img1[zero100],
                   zax=zax[zero100], 
                   zcor=zcor[zero100],
                   zsag=zsag[zero100],
                   izax=izax[zero100], 
                   izcor=izcor[zero100],
                   izsag=izsag[zero100],
                   mean=mom[[1]][zero100],
                   sd=mom[[2]][zero100],                   
                   skew=mom[[3]][zero100],                    
                   mask=mask[zero100],
                   y = truth[zero100])

  # df <- data.frame(val=img1[mask],
  #                  zax=zax[mask], 
  #                  zcor=zcor[mask],
  #                  zsag=zsag[mask],
  #                  izax=izax[mask], 
  #                  izcor=izcor[mask],
  #                  izsag=izsag[mask],
  #                  y = truth[mask])
  # df <- data.frame(val=img1,
  #                  zax=c(zax), 
  #                  zcor=c(zcor),
  #                  zsag=c(zsag),
  #                  izax=c(izax), 
  #                  izcor=c(izcor),
  #                  izsag=c(izsag),
  #                  mask=c(mask),
  #                  y = c(truth))  
  # df <- data.frame(val=c(img1),
  #                  izax=c(izax), 
  #                  izcor=c(izcor),
  #                  izsag=c(izsag),
  #                  y = c(truth))   
  df$over40 <- df$val >= 40
  
  # over0 <- img1 > 0
  # system.time(smooth <- GaussSmoothArray(img1, mask=over0))
  # system.time(smooth <- GaussSmoothArray(img1, mask=mask))
  system.time(smooth <- GaussSmoothArray(img1, mask=zero100))
  # df$smooth <- smooth[mask]
  # df$smooth <- c(smooth)
  df$smooth <- smooth[zero100]

  perc <- 0.1
  N <- nrow(df)
  samp <- sample(1:N, floor(N*perc), replace=FALSE)
  
  train <- df[samp,]
  test <- df[ !(1:N %in% samp), ]
  print(sum(train$y))
  mod <- glm(y ~ . - mask, data=train, family=binomial)
  summary(mod)

  mod2 <- glm(y ~ . - mask - zcor -zsag, 
    data=train, family=binomial)

  mod3 <- glm(y ~ . - mask - smooth, 
    data=train, family=binomial)

  mod4 <- glm(y ~ . - mask - smooth - zcor - zsag, 
    data=train, family=binomial)

  mod5 <- glm(y ~ . - mask - mean - sd - skew, 
    data=train, family=binomial)

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
  mods <- list(mod, mod2, mod3, mod4, mod5)
  mods <- lapply(mods, function(x) {
    x$data <- NULL
    x
  })
  preds <- lapply(mods, function(model) 
    predict(model, newdata=test, type="response") )
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
  outdir <- file.path(dirname(mat$img[irow]), "prob_maps")  
  fname <- basename(mat$img[irow])
  pdf(file.path(outdir, paste0(fname, ".pdf")))
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
    text(1, seq(.8, 0.2, length=length(vars)), vars)
  }
  dev.off()
  
  ### save results to rda  - clear data to free up space
  mod$data <- NULL
  pimg <- img1
  # pimg@.Data[!mask] <- 0
  all.pred <- predict(mod, newdata=df, type="response")
  all.pred <- array(all.pred, dim=dim(pimg@.Data))
  # pimg@.Data[mask] <- all.pred
  pimg@.Data <- all.pred
  rda <- paste0(mat$img[irow], ".rda")
  save(list = c("mods", 
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
  
# }