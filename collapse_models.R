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
basedir <- "/Volumes/DATA_LOCAL/Image_Processing/Test_5"
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

nfiles <- nrow(mat)
irow <- 1
rda <- paste0(mat$img[irow], ".rda")
load(file=rda)
df$fname <- mat$img[irow]
all.df <- df
rm(list="df")
gc()
for (irow in 2:nfiles){
  rda <- paste0(mat$img[irow], ".rda")
  load(file=rda) 
  df$fname <- mat$img[irow]
  all.df <- rbind(all.df, df)
  rm(list="df")
  for (i in 1:10) gc()
  print(irow) 
}


df <- all.df
all.df <- NULL
  perc <- 0.1
  indices <- 1:nrow(df)
  scan.to.drop <- nfiles
  dropscan <- mat$img[scan.to.drop]

  drop.ind <- which(df$fname %in% dropscan )
  indices <- indices[ -drop.ind ]

  N <- length(indices)
  samp <- sample(indices, floor(N*perc), replace=FALSE)
  
  train <- df[samp,]
  test <- df[ !(1:nrow(df) %in% samp), ]

  mod <- glm(y ~ . - fname - mask, data=train, family=binomial)
  summary(mod)
  
  mod$data <- NULL
  rda <- file.path(basedir, "Collaped_Model.rda")
  save(list = c("mod", 
            "scan.to.drop", 
            "dropscan"),
            file=rda)

  #test$pred <- predict(mod, newdata=test, type="response")
  #sum(train$y)
  
  
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

  tt <- NULL
  #### predictions are made!
  test$pred <- pred
  pred <- NULL

  ### data set where the scans were fit on
  within.test <- test[ !test$fname %in% dropscan, ]
  ### getting ROC Curve 
  rpred <- prediction(predictions=within.test$pred, 
    labels=within.test$y)
  rperf <- performance(rpred, "tpr", "fpr")

  outdir <- file.path(basedir, "prob_maps")  
  fname <- "All_Images_NoDropScan"
  pdf(file.path(outdir, paste0(fname, ".pdf")))
    plot(rperf, main=fname)
    abline(v=0.05, col="red")
    abline(v=0.1, col="blue")    
    abline(h=0.9)        
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
    iscan <- 1
    lwd <- ifelse(iscan == scan.to.drop, 3, 1)

    plot(rperfs[[iscan]], col=1, lwd=lwd)
    for (iscan in 1:length(rperfs)){
      lwd <- ifelse(iscan == scan.to.drop, 3, 1)
      plot(rperfs[[iscan]], add=TRUE, col=iscan, lwd=lwd)
    }
    abline(v=0.05, col="red")
    abline(v=0.1, col="blue")
    abline(h=0.9)    
  dev.off()

