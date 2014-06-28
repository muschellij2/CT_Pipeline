library(fslr)
library(cttools)
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
  # source(file.path(progdir, "convert_DICOM.R"))
  
  basedir <<- file.path(homedir, id)
  
}
id <- "205-517"
setup(id)

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

### This is to threshold the values again
rethresh <- FALSE

iscen = as.integer(Sys.getenv("SGE_TASK_ID"))
if (is.na(iscen)) stop("no scenario denoted")

if (rethresh){
  # for (iscen in 1:nrow(scenarios)){
    mask <- scenarios$mask[iscen]
    lower <- scenarios$cutoffs[iscen]
    outdir <- file.path(basedir, lower)
    if (!file.exists(outdir)) {
      system(sprintf('mkdir -p "%s"', outdir))
    }
    for (ifname in fnames){
      stub = nii.stub(basename(ifname))
      addstub = paste0("L", lower)
      if (mask) addstub <- paste0(addstub, "_mask")
      outfile =  file.path(outdir, paste0(stub, "_", addstub))

      fslthresh(file=ifname, outfile=outfile, thresh=lower, 
                mask=mask)
      cat(paste0(ifname, "\n"))
    }
    print(lower)
  # }
}



mask <- scenarios$mask[iscen]
lower <- scenarios$cutoffs[iscen]
outdir <- file.path(basedir, lower)
addstub = paste0("L", lower)
if (mask) addstub <- paste0(addstub, "_mask")

fnames <- list.files(path=outdir, pattern=paste0(addstub, ".nii.gz"), 
                     recursive=FALSE, full.names = TRUE)
fnames <- fnames[!grepl("rigid", fnames)]

ifname <- fnames[1]
ref.img <- fnames[1]
  
  for (iname in 2:length(fnames)){
    ifname <- fnames[iname]

    stub <- paste0("rigid_", basename(ifname))
    ## strip off .nii.gz or .nii
    
    outdir <- dirname(ifname)
    outfile <- file.path(outdir, paste0(stub))
    outmat = paste0(nii.stub(outfile), ".mat")

    ret = flirt(infile=ifname, reffile=ref.img, outfile=outfile,
      omat = outmat,
      dof = 6, intern=FALSE, retimg = FALSE, 
      opts = "-cost leastsq -v" )
    ## used mutualinfo before
    # flirt.wrap(image=ifname, rigid=TRUE, run=TRUE, ref=ref.img, 
    #            mask=mask)
    cat("\n\n")
#     bn <- basename(ifname)
#     fslthresh(image=ifname, autoname=TRUE, outdir=outdir, lower=icut)
  }
  print(icut)
# }
