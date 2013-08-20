library(oro.nifti)
library(animation)
rm(list=ls())
basedir <- "~/CT_Registration/ICES/205-509"
basedir <- path.expand(basedir)

setwd(basedir)

fnames <- list.files(path=basedir, pattern=".nii.gz", 
                     recursive=FALSE, full.names = TRUE)
ani.options(outdir = basedir)
cutoffs <- seq(100, 1000, by=100)


ifname <- fnames[1]
for (ifname in fnames){
  stub <- basename(ifname)
  stub <- gsub("\\.nii\\.gz", "", stub)
  x <- readNIfTI(ifname, 
                 reorient=FALSE)
  
  saveGIF({
  for (icut in cutoffs){
  #   fname <- sprintf("Cutoff_%04.0f.png", icut)
  #   png(fname)
    xx <- x
    xx[xx < icut] <- 0
    xx[xx > 0] <- 1
    orthographic(xx, xyz=c(250, 250, 20), text=icut)
  #   dev.off()
    print(icut)
  }
  }, movie.name= paste0(stub, ".gif"))
  # im.convert("Cutoff_*.png", output="Cutoffs.gif")
  # system('rm Cutoff_*.png')
}