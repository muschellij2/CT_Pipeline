################################
# Written 2013Oct24
# Author: John Muschelli
# Purpose: Compare coefficient estimates from 5 different
# scans, each a different subject in ICES
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

N <- nrow(mat)
irow <- 1
rda <- paste0(mat$img[irow], ".rda")
load(file=rda)
coefs <- coef(mod)
mods <- matrix(NA, ncol=N, nrow=length(coefs))
mods[, 1] <- coefs

for (irow in 2:N){
  rda <- paste0(mat$img[irow], ".rda")
  load(file=rda) 
  coefs <- coef(mod)
  mods[, irow] <- coefs  
}

rownames(mods) <- names(coefs)
colnames(mods) <- basename(mat$img)

outdir <- file.path(basedir, "prob_maps")  

pdf(file.path(outdir, "Coefficient_Image.pdf"))
image(t(mods), yaxt='n')
l <- length(coefs) 
L <- 1+2/l
at <- seq(0, 1, by=L/(l+1))
# at <- at[-length(at)]
axis(2, at=at, labels = FALSE)
text(y = at, par("usr")[1], labels = names(coefs), 
    srt = 0, pos = 2, xpd = TRUE, cex=0.75)
dev.off()