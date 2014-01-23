################################
# Written 2013Oct29
# Author: John Muschelli
# Purpose: Make NIfTI files from DICOMS
# Output: NIfTI images for images and ROIs, and slicethickness image
# Slice thickness image can be used for volume weighting
################################

rm(list=ls())
library(oro.dicom)
library(oro.nifti)
basedir <- "/Volumes/DATA_LOCAL/Image_Processing/Test_Registration"
if (Sys.info()[["user"]] %in% "jmuschel") {
  basedir <- "/dexter/disk2/smart/stroke_ct/ident/Test_Registration"
}
dirs <- list.dirs(basedir, recursive=FALSE, full.names=FALSE)
ptpat <- "\\d\\d\\d-(\\d|)\\d\\d\\d"
ids <- grep(paste0(ptpat, "$"), dirs, value=TRUE)

idir <- 2
## Should we re-run the images (even if png is made?)
redo = TRUE

whichdir = "reoriented"
# whichdir = "FLIRT"
redir = file.path(basedir, whichdir)
imgdir <- file.path(redir, "results")
dir.create(imgdir, showWarnings=FALSE)

fnames <- list.files(path=redir, pattern="^(bws|w|affine).*.nii", 
                     recursive=FALSE, full.names = TRUE)

ifname <- 1

run = data.frame(fname=fnames, stringsAsFactors=FALSE)
run$run = FALSE

if (redo) system(sprintf('rm "%s"/*.png "%s"/*.pdf', imgdir, imgdir))

for (ifname in seq_along(fnames)){
  
  fname <- fnames[ifname]
  imgname <- gsub("\\.nii(|\\.gz)$", "", fname)
  bn <- basename(imgname)
  #bn = gsub("^(w|bws|affine(9|12)_)", "", bn)

  opngname <- pngname <- file.path(imgdir, 
    paste0(bn, ".png"))
  if (!file.exists(pngname) | redo) {
    y <- try(x <- readNIfTI(fname, 
                            reorient=TRUE))
    if (inherits(y, "try-error")) {
      print(paste0(bn, " did not read in"))
      next;
    }
    run$run[ifname] = TRUE
    dx <- dim(x)
    x[x < 0] <- 0
    x[x > 100] <- 100

    x@cal_max <- max(x, na.rm=TRUE)
    x@cal_min <- min(x, na.rm=TRUE)
    
    png(pngname)
    if (all(dim(x) > 1)) {
      orthographic(x, xyz=ceiling(dx/2), 
        text=bn, text.cex=0.8)
    } else {
      plot(0, 0, type='n')
      text(0, 0, bn)
    }
    dev.off()
    # if (opngname != pngname){
    # system(sprintf('mv "%s" "%s"', pngname, opngname))
    # }
  }
  print(ifname)
}

  pdfname <- file.path(imgdir, paste0(whichdir, "_All_Images.pdf"))
  system(sprintf('convert "%s"/*.png "%s"', imgdir, pdfname))

if (whichdir %in% c("FLIRT")){
  pdfname <- file.path(imgdir, 
    paste0(whichdir, "_All_Affine9_Images.pdf"))
  system(sprintf('convert "%s"/affine9*.png "%s"', imgdir, 
    pdfname))  

  pdfname <- file.path(imgdir, 
    paste0(whichdir, "_All_Affine12_Images.pdf"))
  system(sprintf('convert "%s"/affine12*.png "%s"', imgdir, 
    pdfname))    
}