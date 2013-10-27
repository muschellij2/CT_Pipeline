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
  resdir <<- file.path(homedir, "results")
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

# fnames <- list.files(path=basedir, pattern=".nii.gz", 
#                      recursive=FALSE, full.names = TRUE)
# fnames <- fnames[grepl("_CT_", fnames)]
# fnames <- fnames[!grepl("rigid", fnames)]
# 
# ifname <- fnames[1]
masks <- c(FALSE, TRUE)
scenarios <- expand.grid(cutoffs=cutoffs, mask=masks)
iscen <- 1

rethresh <- FALSE

### get all the matrices
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

### make a data set witht he 3x3 matrix as 9 numbers
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

### data is the filename, cutoff, and mask indicator, and 9 numbers
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

pdat <- as.matrix(trans.sd[, cols])
rr <- range(pdat[is.finite(pdat)], na.rm=TRUE)
brks <- seq(from=rr[1], to=rr[2], length=101)
image(x=1:9, y = 1:nrow(pdat), z=t(pdat), xaxt="n", 
      col=heat.colors(100), breaks=brks)
axis(side=3, at=1:9, labels=TRUE)





#revert back to a matrix for tr(Σ^{-1}Σ - I)
matmake <- function(x, n = 3, symm=FALSE) {
  if (!symm) {
    return(matrix(x, nrow=n, ncol=n))
  } else {
    mat <- matrix(NA, nrow=n, ncol=n)
    mat[lower.tri(mat)] <- x
    #mat <- mat + t(mat)
    diag(mat) <- NA
    return(mat)
  }
}

comb <- expand.grid(file=unique(data$file), mask = masks)
comb <- comb[order(comb$file, comb$mask),]
icomb <- 1
## combination of 100-1000
ncuts <- length(cutoffs)
combs <- combn(1:ncuts, 2)
addmat <- matrix(NA, nrow=nrow(comb), ncol = ncol(combs))
cnames <- paste0("V", apply(combs, 2, paste0, collapse="_"))
colnames(addmat) <- cnames
comb <- cbind(comb, addmat)

for (icomb in 1:nrow(comb)){
  imask <- comb$mask[icomb]
  ifile <- comb$file[icomb]
  dat <- data[data$file == ifile & data$mask == imask, ]
  
  
  mat <- as.matrix(dat[,cols])
  ## making them back to matrices
  mats <- alply(mat, 1, matmake)
  scen <- data.frame(t(combs))
  ident <- diag(3)
  scen[, 3] <- NA
  for (iscen in 1:nrow(scen)){
    
    m1 <- mats[[scen[iscen,1]]]
    m2 <- mats[[scen[iscen,2]]]
    mat.diff <- (solve(m1) %*% m2) - ident
        
    scen[iscen,3] <- sum(diag(mat.diff))
  }
  tt <- t(scen[,3])
  comb[icomb, cnames] <- tt
}

pdf(file.path(resdir, paste0(id, "_Skull_Reg_Image.pdf")))
  par(mfrow=c(1,2))
  image(t(as.matrix(comb[ comb$mask , cnames])))
  image(t(as.matrix(comb[ !comb$mask , cnames])))
dev.off()



dat <- as.matrix(comb[,cnames])
## making them back to matrices
mats <- alply(dat, 1, matmake, n = ncuts, symm=TRUE)
## make the picture of matrix differences
ncomb <- nrow(comb)
stopifnot((ncomb %% 2) == 0)
brks <- seq(-5, 5, by=1)
cols <- heat.colors(length(brks)-1)
primage <- function(M, N, mask, brks, cols, labs, ...) {
  image(x=1:N, y= 1:N, z=t(M[N:1,]), xaxt="n", yaxt="n", 
        main=ifelse(mask, "Mask", "Raw"),
        col=cols, breaks=brks, 
        ylab="Cutoffs", xlab="Cutoffs", ...)  
  axis(side=1, at=1:N, labels=labs, cex.axis=0.5)
  axis(side=2, at=N:1, labels=labs, cex.axis=0.5, las=1)  
}
pdf(file.path(resdir, paste0(id, "_Skull_Reg_Compare_Mask.pdf")),
    height=3.5, width=7)
  for (icomb in seq(from=1, to=ncomb-1, by=2)){
    par(mfrow=c(1,2))
    M1 <- mats[[icomb]]
    M2 <- mats[[icomb+1]]
    fname <- unique(comb[icomb:(icomb+1),"file"])
    fname <- as.character(fname)
    fname <- gsub(paste0(id, "_"), "", fname)
    N <- nrow(M1)

    primage(M1, N=N, mask=comb$mask[icomb], brks, cols, labs=cutoffs)
    primage(M2, N=N, mask=comb$mask[icomb+1], brks, cols, labs=cutoffs)
    
    title(outer=TRUE, main=c("", "", fname))
  }
dev.off()
