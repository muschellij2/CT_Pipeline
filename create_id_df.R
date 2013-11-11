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
# imgs <- unlist(mat)
# imgs <- paste0(imgs, ".nii.gz")
# stopifnot(all(file.exists(imgs)))

# irow <- grep("102323", mat$img)

irow <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if (is.na(irow)) irow <- 4

### test to see if all files exist
check <- mat[irow,]
check <- paste0(check, ".nii.gz")
stopifnot(all(file.exists(check)))

# for (irow in 1:nrow(mat)){
  ## manual segmentation
  truth <- readNIfTI(mat$roi[irow], reorient=FALSE)
  
  true.mask <- truth > 0
  stopifnot(all(truth %in% c(0, 1)))

  img1 <- readNIfTI(mat$img[irow], reorient=FALSE)

  vals <- img1[true.mask]
  vals <- data.frame(val=vals)

  outdir <- file.path(basedir, "hist_clot_vals") 
  fname <- basename(mat$img[irow])  
  pdf(file.path(outdir, paste0(fname, ".pdf")))
    g <- ggplot(vals, aes(x=val)) + 
      geom_histogram(aes(y=..density..), fill="red") + 
      xlab("Value") + ylab("Density") + 
      ggtitle(paste0("Distribution of values for Hemorrhage\n", fname)) 
    g2 <- g + geom_line(stat="density")
    print(g)
    print(g2)
  dev.off()
  # system.time(zs <- zscore2(img1))
  system.time(izax <- zscore(img1))
  system.time(izcor <- zscore(img1, margin=2))
  system.time(izsag <- zscore(img1, margin=1))


  ss.mask <- readNIfTI(mat$ss.img[irow], reorient=FALSE)
  mask <- ss.mask > 0

#  zero100 <- img1 >= 0 & img1 <= 100
### masking out non-brain 
  ss.img <- img1
  ss.img[!mask] <- 0

  system.time(mom <- image_moment(img1, mask=mask, 
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
  


  # stopifnot(any(truth[!zero100] == 1))


  # df <- data.frame(val=img1[zero100],
  #                  zax=zax[zero100], 
  #                  zcor=zcor[zero100],
  #                  zsag=zsag[zero100],
  #                  izax=izax[zero100], 
  #                  izcor=izcor[zero100],
  #                  izsag=izsag[zero100],
  #                  mean=mom[[1]][zero100],
  #                  sd=mom[[2]][zero100],                   
  #                  skew=mom[[3]][zero100],                    
  #                  mask=mask[zero100],
  #                  y = truth[zero100])

  df <- data.frame(val=img1[mask],
                   zax=zax[mask], 
                   zcor=zcor[mask],
                   zsag=zsag[mask],
                   izax=izax[mask], 
                   izcor=izcor[mask],
                   izsag=izsag[mask],
                   mean=mom[[1]][mask],
                   sd=mom[[2]][mask],                   
                   skew=mom[[3]][mask],                    
                   y = truth[mask])
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
  pixdims <- img1@pixdim[2:4]
  # system.time(smooth <- GaussSmoothArray(img1, 
  #   voxdim=pixdims,
  #   mask=mask))
  # df$smooth <- smooth[mask]

  img1[!mask] <- 0
  sigma <- 1
  for (sigma in c(1, 3, 5, 10, 20)) {
    fname <- basename(mat$img[irow])
    fname <- gsub(".nii.gz", "", fname, fixed=TRUE)
    outfile <- file.path(basedir, "Smoothed_Images", fname)
    outfile <- sprintf('%s_%s', outfile, sigma)
    fslsmooth(file=mat$img[irow], 
      mask = mat$ss.img[irow], 
      outfile = outfile,
      sigma=sigma, 
      local=FALSE)
    simg <- readNIfTI(outfile, reorient=FALSE)
    simg[!mask] <- 0
    df[, paste0("smooth", sigma)] <- simg[mask]
  }

  # df$smooth <- smooth[mask]
  # df$smooth <- c(smooth)
  rda <- paste0(mat$img[irow], "_Data.rda")
    save(list = "df",
            file=rda)
