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
subset <- sample(1:ncol(combos), 1)
ind <- 1:nfiles
train.ind <- ind %in% combos[, subset]
test.ind <- !train.ind

  rda <- paste0(mat$img[1], "_Models.rda")
  load(file=rda)
  run.mods <- mods

rda <- file.path(basedir, "Collaped_Data.rda")
load(file=rda)

  perc <- 0.1
  indices <- 1:nrow(df)
  dropscan <- mat$img[test.ind]

  drop.ind <- which(df$fname %in% dropscan )
  indices <- indices[ -drop.ind ]

  N <- length(indices)
  samp <- sample(indices, floor(N*perc), replace=FALSE)
  
  train <- df[samp,]
  test <- df[ !(1:nrow(df) %in% samp), ]


  col.mods <- lapply(run.mods, function(mod){
    cmod <- glm(formula=mod$formula, data=train, family=binomial)
    cmod <- scrape.mod(cmod)
    cmod
    })

  # mod <- glm(y ~ . - fname, data=train, family=binomial, model=FALSE)
  # summary(mod)
  # mod <- scrape.mod(mod)

  # mod2 <- glm(y ~ . -fname - zcor -zsag, 
  #   data=train, family=binomial)
  # mod2 <- scrape.mod(mod2)

  # mod3 <- glm(y ~ . -fname - smooth, 
  #   data=train, family=binomial)
  # mod3 <- scrape.mod(mod3)

  # mod4 <- glm(y ~ . -fname - smooth - zcor - zsag, 
  #   data=train, family=binomial)
  # mod4 <- scrape.mod(mod4)

  # mod5 <- glm(y ~ . -fname - mean - sd - skew, 
  #   data=train, family=binomial)
  # mod5 <- scrape.mod(mod5)

  # mods <- list(mod, mod2, mod3, mod4, mod5)

  # mod2 <- mod3 <- mod4 <- mod5 <- NULL
  mods <- lapply(col.mods, scrape.mod)
  rda <- file.path(basedir, "Collaped_Model.rda")
  save(list = c("mods", 
            "train.ind",
            "samp", 
            "dropscan"),
            file=rda)

  #test$pred <- predict(mod, newdata=test, type="response")
  #sum(train$y)
  
  
  ##########################
  ###Getting predictions####
  ##Faster than predict#####
  ##########################
  # get.stuff <- function(mod){
  #   var.used <- attr(mod$terms, "term.labels")
  #   var.classes <- attr(mod$terms, "dataClasses")
  #   var.classes <- var.classes[var.used]
  # }
  # var.classes <- get.stuff(mod)


  # ### checking to see if we used factors or not in the model
  # stopifnot(all(var.classes %in% c('numeric', 'logical')))
  # coefs <- coef(mod)
  # have.log <- var.classes %in% "logical"
  # if (any(have.log)){
  #   vars <- names(var.classes)[have.log]
  #   for (ivar in vars){
  #     names(coefs) <- gsub(paste0(ivar, "TRUE"), ivar, names(coefs))
  #   }
  # }
  # ## just getting the intercept - not needed in matrix mult
  # intercept <- ifelse("(Intercept)" %in% names(coefs), 
  #   coefs["(Intercept)"], 
  #   0)
  # coefs <- coefs[ !names(coefs) %in% "(Intercept)" ]
  # tt <- test[, names(coefs)]

  # ### getting predictions - prob slower than matrix multiplication
  # ### but more control
  # for (icoef in names(coefs)){
  #   tt[, icoef] <- coefs[names(coefs) == icoef] * tt[, icoef]
  #   print(icoef)
  # }
  # pred <- rowSums(tt)
  # pred <- pred + intercept
  # pred <- exp(pred)/(1+exp(pred))

################################################################
  # # tt <- NULL
  # #### predictions are made!
  # preds <- sapply(mods, pred.prob, test=test )
  # nmods <- ncol(preds)
  # cn <- colnames(preds) <- paste0("mod", 1:nmods)
  # test <- cbind(test, preds)

  # test <- test[, c("fname", "y", cn), with=FALSE]
  # test <- as.data.frame(test, stringsAsFactors=FALSE)
  # # test$pred <- pred.prob(mod, test)
  # # pred <- NULL

  # # ### data set where the scans were fit on
  # # within.test <- test[ !test$fname %in% dropscan, ]
  # # ### getting ROC Curve 
  # # rpred <- prediction(predictions=within.test$pred, 
  # #   labels=within.test$y)
  # # rperf <- performance(rpred, "tpr", "fpr")

  # # fname <- "All_Images_NoDropScan"
  # # pdf(file.path(outdir, paste0(fname, ".pdf")))
  # #   plot(rperf, main=fname)
  # #   abline(v=0.05, col="red")
  # #   abline(v=0.1, col="blue")    
  # #   abline(h=0.9)        
  # # dev.off()

  # #### getting scan-specific ROC Curve
  # # dropscan 
  # outdir <- file.path(basedir, "prob_maps")  

  # oos <- test[ test$fname %in% dropscan, ]

  # fname <- "Out_of_Sample_ROC"
  # pdf(file.path(outdir, paste0(fname, ".pdf")))
  # imod <- 1
  # for (imod in 1:length(mods)){
  #   modcol <- paste0("mod", imod)

  #   ndrop <- length(dropscan)

  #   rpreds <- vector("list", length=ndrop)
  #   for (ifname in 1:ndrop){
  #     tt <- oos[ oos$fname %in% dropscan[ifname], ]
  #     rpreds[[ifname]] <- prediction(predictions=tt[, modcol], 
  #       labels=tt$y)
  #     print(ifname)
  #   }
  #   rperfs <- lapply(rpreds, performance, "tpr", "fpr")

  #   iscan <- 1
  #   # lwd <- ifelse(iscan == scan.to.drop, 3, 1)
  #   lwd <- 1
  #   # par(mfrow=c(2, 1))
  #   par(mfrow=c(1, 2))
  #   plot(rperfs[[iscan]], col=1, lwd=lwd)
  #   for (iscan in 1:length(rperfs)){
  #     plot(rperfs[[iscan]], add=TRUE, col=iscan)
  #   }
  #   abline(v=0.05, col="red")
  #   abline(v=0.1, col="blue")
  #   abline(h=0.9)


  #   plot(rperfs[[iscan]], col=1, lwd=lwd, xlim=c(0, 0.1))
  #   for (iscan in 1:length(rperfs)){
  #     plot(rperfs[[iscan]], add=TRUE, col=iscan, xlim=c(0, 0.1))
  #   }
  #   abline(v=0.05, col="red")
  #   abline(v=0.1, col="blue")
  #   abline(h=0.9)
  # }
  
  # dev.off()

  # ### in sample 
  # insamp <- test[ !test$fname %in% dropscan, ]

  # fname <- "In_Sample_ROC"
  # pdf(file.path(outdir, paste0(fname, ".pdf")))
  # imod <- 1
  # for (imod in 1:length(mods)){
  #   modcol <- paste0("mod", imod)

  #   ndrop <- length(dropscan)

  #   rpreds <- vector("list", length=ndrop)
  #   for (ifname in 1:ndrop){
  #     tt <- insamp[ insamp$fname %in% dropscan[ifname], ]
  #     rpreds[[ifname]] <- prediction(predictions=tt[, modcol], 
  #       labels=tt$y)
  #     print(ifname)
  #   }
  #   rperfs <- lapply(rpreds, performance, "tpr", "fpr")

  #   iscan <- 1
  #   # lwd <- ifelse(iscan == scan.to.drop, 3, 1)
  #   lwd <- 1
  #   # par(mfrow=c(2, 1))
  #   par(mfrow=c(1, 2))
  #   plot(rperfs[[iscan]], col=1, lwd=lwd)
  #   for (iscan in 1:length(rperfs)){
  #     plot(rperfs[[iscan]], add=TRUE, col=iscan)
  #   }
  #   abline(v=0.05, col="red")
  #   abline(v=0.1, col="blue")
  #   abline(h=0.9)


  #   plot(rperfs[[iscan]], col=1, lwd=lwd, xlim=c(0, 0.1))
  #   for (iscan in 1:length(rperfs)){
  #     plot(rperfs[[iscan]], add=TRUE, col=iscan, xlim=c(0, 0.1))
  #   }
  #   abline(v=0.05, col="red")
  #   abline(v=0.1, col="blue")
  #   abline(h=0.9)
  # }
  
  # dev.off()  
################################################################



# rpreds <- vector("list", length=ndrop)
#   for (ifname in 1:ndrop){
#     tt <- test[ test$fname %in% mat$img[ifname], ]
#     rpreds[[ifname]] <- prediction(predictions=tt$pred, labels=tt$y)
#     print(ifname)
#   }
#   rperfs <- lapply(rpreds, performance, "tpr", "fpr")

#   iscan <- 1
#   lwd <- ifelse(iscan == scan.to.drop, 3, 1)

#   plot(rperfs[[iscan]], col=1, lwd=lwd)
#   for (iscan in 1:length(rperfs)){
#     lwd <- ifelse(iscan == scan.to.drop, 3, 1)
#     plot(rperfs[[iscan]], add=TRUE, col=iscan, lwd=lwd)
#   }
#   abline(v=0.05, col="red")
#   abline(v=0.1, col="blue")
#   abline(h=0.9)    
#   dev.off()