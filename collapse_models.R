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

N <- nrow(mat)
irow <- 1
rda <- paste0(mat$img[irow], ".rda")
load(file=rda)
df$fname <- mat$img[irow]
all.df <- df

for (irow in 2:N){
  rda <- paste0(mat$img[irow], ".rda")
  load(file=rda) 
  df$fname <- mat$img[irow]
  all.df <- rbind(all.df, df)
}

df <- all.df
all.df <- NULL
  perc <- 0.1
  N <- nrow(df)
  samp <- sample(1:N, floor(N*perc), replace=FALSE)
  
  train <- df[samp,]
  test <- df[ !(1:N %in% samp), ]

  mod <- glm(y ~ . - fname, data=train, family=binomial)
  summary(mod)
  
  test$pred <- predict(mod, newdata=test, type="response")
  sum(train$y)
  
  
  ##########################
  ###Getting predictions####
  ##Faster than predict#####
  ##########################
  get.stuff <- function(mod){
    var.used <- attr(mod$terms, "term.labels")
    var.classes <- attr(mod$terms, "dataClasses")
    var.classes <- var.classes[var.used]
  }
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
  ## just getting the intercept - not needed in matrix mult
  intercept <- ifelse("(Intercept)" %in% names(coefs), 
    coefs["(Intercept)"], 
    0)
  coefs <- coefs[ !names(coefs) %in% "(Intercept)" ]
  tt <- test[, names(coefs)]

  ### getting predictions - prob slower than matrix multiplication
  ### but more control
  for (icoef in names(coefs)){
    tt[, icoef] <- coefs[names(coefs) == icoef] * tt[, icoef]
    print(icoef)
  }
  pred <- rowSums(tt)
  pred <- pred + intercept
  pred <- exp(pred)/(1+exp(pred))

  #### predictions are made!
  test$pred <- pred


  ### getting ROC Curve 
  rpred <- prediction(predictions=pred, labels=test$y)
  rperf <- performance(rpred, "tpr", "fpr")

  outdir <- file.path(dirname(mat$img[irow]), "prob_maps")  
  fname <- "All_Images"
  pdf(file.path(outdir, paste0(fname, ".pdf")))
    plot(rperf, main=fname)
  dev.off()

  #### getting scan-specific ROC Curve
  nfiles <- nrow(mat)
  rpreds <- vector("list", length=nfiles)
  for (ifname in 1:nfiles){
    tt <- test[ test$fname %in% mat$img[ifname], ]
    rpreds[[ifname]] <- prediction(predictions=tt$pred, labels=tt$y)
    print(ifname)
  }
  rperfs <- lapply(rpreds, performance, "tpr", "fpr")

  fname <- "All_Images_PerScan_ROC"
  pdf(file.path(outdir, paste0(fname, ".pdf")))
    plot(rperfs[[1]], col=1)
    for (iscan in 1:length(rperfs)){
      plot(rperfs[[iscan]], add=TRUE, col=iscan)
    }

  dev.off()

