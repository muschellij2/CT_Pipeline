library(oro.nifti)
rm(list=ls())
setup <- function(id){
  username <- Sys.info()["user"][[1]]
  
  cluster=FALSE
  if (username == "muschellij2"){
    # rootdir <- "/Volumes/DATA/New_Age_Test"
    rootdir <- "~/CT_Registration"
  } else {
    rootdir <- "/dexter/disk2/smart/stroke_ct/ident"
    cluster =TRUE;
  }
  rootdir <- path.expand(rootdir)
  
  ss <- as.numeric(strsplit(id, "-")[[1]][2])
  if (ss > 4000){
    study <- "CLEAR_III"
    dpath <- file.path("CLEAR", "CLEAR III")
  } else if (ss > 300 & ss < 500){
    dpath <- study <- "MISTIE"
  } else if (ss > 500 & ss < 4000) {
    dpath <- study <- "ICES" 
  }
  
  
  rootdir <<- path.expand(rootdir)
  homedir <<- file.path(rootdir, study)
  homedir <<- path.expand(homedir)
  
  #progdir <- file.path(dirname(basedir), "programs")
  progdir <<- file.path(rootdir, "programs")
  source(file.path(progdir, "convert_DICOM.R"))
  
  basedir <<- file.path(homedir, id)
  
}
id <- "205-509"
setup(id)
source(file.path(progdir, "fslthresh.R"))

setwd(basedir)


cutoffs <- seq(100, 1000, by=100)
# cutoffs <- 300
iname <- 1

fnames <- list.files(path=basedir, pattern=".nii.gz", 
                     recursive=FALSE, full.names = TRUE)
fnames <- fnames[grepl("_CT_", fnames)]
fnames <- fnames[!grepl("rigid", fnames)]

ifname <- fnames[1]
masks <- c(FALSE, TRUE)
scenarios <- expand.grid(cutoffs=cutoffs, mask=masks)
iscen <- 1

rethresh <- FALSE

files <- list.files(path=basedir, 
                     pattern="rigid.*.mat", 
                     recursive=TRUE, full.names = TRUE)
fnames <- basename(files)
fnames <- gsub(".mat", "", fnames)

ss <- strsplit(fnames, split="_")
file.base <- t(sapply(ss, function(x){
  x <- x[-1]
  lx <- length(x)
  havemask <- any(grepl( "mask", x[(lx-1):lx]))
  if (havemask) x <- x[-length(x)]
  lx <- length(x)
  xname <- paste0(x[1:(lx-1)], sep="", collapse="_")
  return(c(xname, x[length(x)]))
}))

data <- data.frame(file=file.base[,1], stringsAsFactors=TRUE)
data$cut <- as.numeric(gsub("L", "", file.base[,2]))
data$fname <- paste0(fnames, ".mat")
data$ID <- as.numeric(data$file)
data$mask <- grepl("_mask", data$fname)
data$fullfile <- files
data <- data[order(data$ID, data$cut, data$mask),]



### readin in matrices and getting averages
dd <- sapply(data$fullfile, function(x) {
  tab <- as.matrix(read.table(x))
  tab <- tab[1:3, 1:3]
  tab <- c(tab)
  return(tab)
})
dd <- t(dd)
rownames(dd) <- NULL
colnames(dd) <- sapply(1:3, function(x) paste0("V", 1:3, x))

data <- cbind(data, dd)
cols <- c(sapply(1:3, function(x) paste0("V", 1:3, x)))


trans.sd <- aggregate(cbind(V11, V21, V31, V12, V22, V32, V13, V23, V33) ~ ID + mask, 
  data=data, sd)
trans.mn <- aggregate(cbind(V11, V21, V31, V12, V22, V32, V13, V23, V33) ~ ID + mask, 
  data=data, mean)
### getting coeff of variation
trans.cv <- trans.mn[, cols] / trans.sd[, cols]
trans.cv <- cbind(trans.sd[, 1:2], trans.cv)

### round these to look at
# trans.sd[, cols ] <- round(trans.sd[, cols ], 3)
# trans.mn[, cols ] <- round(trans.mn[, cols ], 3)

trans.sd <- trans.sd[order(trans.sd$ID, trans.sd$mask), ]
trans.mn <- trans.mn[order(trans.mn$ID, trans.mn$mask), ]
trans.cv <- trans.cv[order(trans.cv$ID, trans.cv$mask), ]



