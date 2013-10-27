library(oro.nifti)
rm(list=ls())

basedir <- "~/CT_Registration/ICES"
basedir <- path.expand(basedir)
setwd(basedir)
ptpat <- "\\d\\d\\d-(\\d|)\\d\\d\\d"

dirs <- list.dirs(basedir, full.names=FALSE, recursive=FALSE)
dirs <- dirs[grepl(paste0(ptpat, '$'), dirs)]
idir <- 2
## Should we re-run the images (even if png is made?)
redo = FALSE
# dirs <- dirs[!grepl("232-514", dirs)]
# dirs <- dirs[!grepl("232-513", dirs)]
# dirs <- dirs[!grepl("232-512", dirs)]
# dirs <- dirs[!grepl("225-524", dirs)]
# dirs <- dirs[!grepl("225-523", dirs)]

for (idir in seq_along(dirs)){
  path <- dirs[idir]
  imgdir <- file.path(path, "Images")
  dir.create(imgdir, showWarnings=FALSE)
  fnames <- list.files(path=path, pattern=".nii.gz", 
                     recursive=FALSE, full.names = TRUE)
  # print(fnames)
  ifname <- 1
  for (ifname in seq_along(fnames)){

    fname <- fnames[ifname]
    imgname <- gsub("\\.nii\\.gz$", "", fname)
    bn <- basename(imgname)
    opngname <- pngname <- file.path(imgdir, gsub("\\.nii\\.gz$", ".png", basename(fname)))
    pngname <- gsub("%", "", pngname)
    if (!file.exists(pngname) | redo) {
      y <- try(x <- readNIfTI(fname, 
                   reorient=FALSE))
      if (class(y) == "try-error") next;
      CT <- grepl("_CT_", fname)
      dx <- dim(x)
      if (CT) {
        x[x < 0] <- 0
        x[x > 100] <- 100
      }
      x@cal_max <- max(x, na.rm=TRUE)
      x@cal_min <- min(x, na.rm=TRUE)

      bn <- gsub(paste0("^", ptpat, "_"), "", bn)
      png(pngname)
        if (all(dim(x) > 1)) {
          orthographic(x, xyz=ceiling(dx/2), text=bn, text.cex=0.8)
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

  id <- gsub(paste0(".*(", ptpat, ")"), "\\1", path)
  stopifnot(all(grepl(ptpat, id)))
  pdfname <- file.path(basedir, "All_Images", paste0(id, "_ortho_images.pdf"))
  system(sprintf('convert "%s"/*.png "%s"', imgdir, pdfname))
    # dev.off()
  print(id)
}