rm(list=ls())
library(oro.nifti)
library(stringr)
library(R.matlab)

username <- Sys.info()["user"][[1]]

basedir <- file.path("/Users", username, 
                     "Dropbox/CTR/DHanley/CT_Registration/ICES")
source(file.path(basedir, "programs", "fslhd.R"))                    

setwd(basedir)
prepost <- read.csv("Matched_Pre_Post_Op_Scans.csv", 
                    stringsAsFactors=FALSE)

pp <- prepost[, c("patientName", "Type", "Scan_to_Take")]
res <- reshape(pp, direction="wide", idvar="patientName", timevar="Type")
res <- res[res$"Scan_to_Take.Pre-Op" != "",]
res <- res[res$"Scan_to_Take.Post-Op" != "",]
rownames(res) <- NULL
iid <- 1

skullstrip <- TRUE
ssadd <- NULL
if (skullstrip) {
  res[, c("Scan_to_Take.Pre-Op", "Scan_to_Take.Post-Op")] <- 
  sapply(res[, c("Scan_to_Take.Pre-Op", "Scan_to_Take.Post-Op")], 
         gsub, pattern=".nii.gz", replacement="_SS.nii.gz", 
         fixed=TRUE)
  ssadd <- "Skull_Stripped"
}
ssadd2 <- ssadd
if (skullstrip) ssadd2 <- paste0(ssadd2, "_")

# ids <- (nrow(res)-1):1
ids <- 1:nrow(res)
res2 <- res


for (iid in ids){
  id <- paste(substr(res$patientName[iid], 1, 3), substr(res$patientName[iid], 4, 6), sep="-")
  iddir <- file.path(basedir, id, "Separated")
  
  regdir <- file.path(iddir, paste("SPM", ssadd, sep=.Platform$file.sep))
  iddir <- file.path(basedir, id, paste("Separated", ssadd, sep=.Platform$file.sep))
  
  if (!file.exists(regdir)) {
    system(sprintf('mkdir "%s"', regdir))
  }  else {
    system(sprintf('rm "%s"/*.nii', regdir))    
  }
  
  pre <- res$"Scan_to_Take.Pre-Op"[iid]
  refscan <- str_trim(file.path(iddir, pre))
  xx <- readNIfTI(refscan, reorient=FALSE)
  if (min(xx) < -1024) {
    xx[xx < -1024] <- -1024
    refscan <- gsub(".nii.gz", "_Added", refscan)
    xx@cal_max <- max(xx)
    xx@cal_min <- min(xx)
    writeNIfTI(xx, file=refscan)
    refscan <- paste0(refscan, ".nii.gz")
    print(paste("Making new refscan", refscan))
    
  }
  
  post <- res$"Scan_to_Take.Post-Op"[iid]
  inscan <- str_trim(file.path(iddir, post))
  xx <- readNIfTI(inscan, reorient=FALSE)
  if (min(xx) < -1024) {
    xx[xx < -1024] <- -1024
    inscan <- gsub(".nii.gz", "_Added", inscan)
    xx@cal_max <- max(xx)
    xx@cal_min <- min(xx)
    writeNIfTI(xx, file=inscan)
    inscan <- paste0(inscan, ".nii.gz")
    print(paste("Making new inscan", inscan))
  }
  
  
  pre <- res2$"Scan_to_Take.Pre-Op"[iid] <-  basename(refscan)
  post <- res2$"Scan_to_Take.Post-Op"[iid] <-  basename(inscan)

  rrefscan <- str_trim(file.path(regdir, pre ))
  rinscan  <- str_trim(file.path(regdir, post))
  
  system(sprintf('cp "%s" "%s"', refscan, regdir))
  system(sprintf('gunzip -f "%s"', rrefscan))
  
  system(sprintf('cp "%s" "%s"', inscan, regdir))
  system(sprintf('gunzip -f "%s"', rinscan))
  
  print(iid)
  
}


res2$patientName <- paste(substr(res2$patientName, 1, 3), 
                          substr(res2$patientName, 4, 6), sep="-")
colnames(res2) <- c("patientName", "PreOp", "PostOp")
res2$path <- sapply(res2$patientName, 
                    function(x) 
                      file.path(basedir, x, "Separated", "SPM"))
res2$PreOp <- gsub(".nii.gz", ".nii", res2$PreOp, fixed=TRUE)
res2$PostOp <- gsub(".nii.gz", ".nii", res2$PostOp, fixed=TRUE)
res2$PreOp <- file.path(res2$path, res2$PreOp)
res2$PostOp <- file.path(res2$path, res2$PostOp)
res2$path <- NULL
res2 <- res2[complete.cases(res2),]
writeMat(file.path(basedir, paste0(ssadd2, "Scans.mat"), res=res2)
write.csv(res2, file=file.path(basedir, paste0(ssadd2, "Pre_Op_Post_Op_Reg_Scans.csv")))
