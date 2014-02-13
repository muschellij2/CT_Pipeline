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

## all combinations
nfiles <- nrow(mat)
combos <- combn(nrow(mat), floor(nrow(mat)/2))
ind <- 1:nfiles

rda <- file.path(basedir, "Collaped_Data.rda")
x <- load(file=rda)

rda <- file.path(basedir, "Collaped_Model.rda")
x2 <- load(file=rda)

test.ind <- !train.ind

  test <- df[ !(1:nrow(df) %in% samp), ]

  #test$pred <- predict(mod, newdata=test, type="response")
  #sum(train$y)


  # tt <- NULL
  #### predictions are made!
  system.time(preds <- sapply(mods, pred.prob, test=test ))
  nmods <- ncol(preds)
  cn <- colnames(preds) <- paste0("mod", 1:nmods)
  preds <- data.table(preds)
  test <- cbind(test, preds)
  preds <- NULL

  test <- test[, c("fname", "y", cn), with=FALSE]
  test <- as.data.frame(test, stringsAsFactors=FALSE)
  # test$pred <- pred.prob(mod, test)
  # pred <- NULL

  #### getting scan-specific ROC Curve
  # dropscan 
  print.rocs <- function(runsamp, fnames, imod){
      modcol <- paste0("mod", imod)

      nfnames <- length(fnames)

      rpreds <- vector("list", length=nfnames)
      for (ifname in 1:nfnames){
        tt <- runsamp[[ fnames[ifname] ]]
        rpreds[[ifname]] <- prediction(predictions=tt[, modcol], 
          labels=tt$y)
        print(ifname)
      }
      rperfs <- lapply(rpreds, performance, "tpr", "fpr")

      iscan <- 1
      # lwd <- ifelse(iscan == scan.to.drop, 3, 1)
      lwd <- 1
      # par(mfrow=c(2, 1))
      par(mfrow=c(1, 2))
      plot(rperfs[[iscan]], col=1, lwd=lwd)
      for (iscan in 1:length(rperfs)){
        plot(rperfs[[iscan]], add=TRUE, col=iscan)
      }
      abline(v=0.05, col="red")
      abline(v=0.1, col="blue")
      abline(h=0.9)


      plot(rperfs[[iscan]], col=1, lwd=lwd, xlim=c(0, 0.1))
      for (iscan in 1:length(rperfs)){
        plot(rperfs[[iscan]], add=TRUE, col=iscan, xlim=c(0, 0.1))
      }
      abline(v=0.05, col="red")
      abline(v=0.1, col="blue")
      abline(h=0.9)
  }

  outdir <- file.path(basedir, "prob_maps")  

  oos <- test[ test$fname %in% dropscan, ]
  oos <- split(oos, oos$fname)

  fname <- "Out_of_Sample_ROC"
  for (imod in 1:length(mods)){
    png(file.path(outdir, paste0(fname, imod, ".png")))
    print.rocs(oos, dropscan, imod)
    dev.off()
  }      

  rm(list="oos")
  ### in sample 
  sampscan <- mat$img[! mat$img %in% dropscan]
  insamp <- test[ !test$fname %in% dropscan, ]
  insamp <- split(insamp, insamp$fname)

  fname <- "In_Sample_ROC"
  for (imod in 1:length(mods)){
    png(file.path(outdir, paste0(fname, imod, ".png")))
    print.rocs(insamp, sampscan, imod)
    dev.off()
  }      

  rm(list="test")

  train <- df[samp,]
  
## Training data

  system.time(preds <- sapply(mods, pred.prob, test=train ))
  nmods <- ncol(preds)
  cn <- colnames(preds) <- paste0("mod", 1:nmods)
  train <- cbind(train, preds)
  preds <- NULL

  train <- train[, c("fname", "y", cn), with=FALSE]
  train <- as.data.frame(train, stringsAsFactors=FALSE)


  insamp <- train[ !train$fname %in% dropscan, ]
  insamp <- split(insamp, insamp$fname)

  fname <- "In_Sample_Train_ROC"
  for (imod in 1:length(mods)){
    png(file.path(outdir, paste0(fname, imod, ".png")))
    print.rocs(insamp, sampscan, imod)
    dev.off()
  }      


  # pdf(file.path(outdir, paste0(fname, ".pdf")))
  #   print.rocs(insamp, sampscan)
  # dev.off()

