rm(list=ls())
library(oro.nifti)
library(stringr)
library(R.matlab)

username <- Sys.info()["user"][[1]]

basedir <- file.path("/Users", username, 
                     "Dropbox/CTR/DHanley/CT_Registration/ICES")
setwd(basedir)
res2 <- read.csv(file=file.path(basedir, 
                "Pre_Op_Post_Op_Reg_Scans.csv"), 
                 stringsAsFactors=FALSE)
res2$X <- NULL
res2$PreOp <-  paste0(res2$PreOp, ".gz")


scan <- paste0("r", basename(res2$PostOp))
res2$PostOp <- file.path(dirname(res2$PostOp), scan)


iid <- 1
for (iid in 1:nrow(res2)){
  cmd <- sprintf('/usr/local/fsl/bin/fslview.app/Contents/MacOS/fslview "%s" "%s"', 
                 res2$PreOp[iid], res2$PostOp[iid])
  system(cmd)
}