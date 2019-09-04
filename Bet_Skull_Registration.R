rm(list=ls())
library(cttools)
library(fslr)
options(matlab.path='/Applications/MATLAB_R2014b.app/bin')

setup <- function(id, study = "Registration"){
  username <- Sys.info()["user"][[1]]

  cluster=FALSE
  if (username %in% c("muschellij2", "johnmuschelli")){
    # rootdir <- "/Volumes/DATA/New_Age_Test"
    rootdir <- "~/CT_Registration"
  } else {
    rootdir <- "/dexter/disk2/smart/stroke_ct/ident"
    cluster =TRUE;
  }
    rootdir <- path.expand(rootdir)

  # ss <- as.numeric(strsplit(id, "-")[[1]][2])
  # if (ss > 4000){
  #   study <- "CLEAR_III"
  #   dpath <- file.path("CLEAR", "CLEAR III")
  # } else if (ss > 300 & ss < 500){
  #   dpath <- study <- "MISTIE"
  # } else if (ss > 500 & ss < 4000) {
  #   dpath <- study <- "ICES" 
  # }


  rootdir <<- path.expand(rootdir)
  homedir <<- file.path(rootdir, study)
  homedir <<- path.expand(homedir)

#progdir <- file.path(dirname(basedir), "programs")
  progdir <- file.path(rootdir, "programs")
  # source(file.path(progdir, "convert_DICOM.R"))
  # source(file.path(progdir, "fslhd.R"))

  basedir <<- file.path(homedir, id)

}


#### setting up if things are on the cluster or not
study = "Registration"

setup(study, study=study)

ids = list.dirs(homedir, recursive=FALSE, full.names=FALSE)
ids = basename(ids)
ids = grep("\\d\\d\\d-(\\d|)\\d\\d\\d", ids, value=TRUE)
length(ids)

### initial setup

iid <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if (is.na(iid)) iid <- 1

id <- ids[iid]
setup(id, study = study)

setwd(basedir)


# cutoffs <- 300
iname <- 1

fnames <- list.files(path=basedir, pattern=".nii.gz", 
                       recursive=FALSE, full.names = TRUE)
fnames <- fnames[grepl("_CT_", fnames)]
fnames <- fnames[!grepl("rigid", fnames)]

intensity = .01
stub = paste0(nii.stub(basename(fnames)), 
  sprintf("_SS_Mask_%0.2g", intensity), ".nii.gz")
ss = file.path(dirname(fnames), "Skull_Stripped", stub)

df = data.frame(img = fnames, 
  ss = ss, 
  stringsAsFactors=FALSE)
exist = apply(df[, c("img", "ss")], 1, file.exists)
exist = apply(exist, 2, all)

df = df[exist, ]
ifname <- df$img[1]
masks <- c(FALSE, TRUE)
scenarios <- expand.grid(mask=masks)
iscen <- 2

### This is to threshold the values again
rethresh <- FALSE

# for (mask in masks){
mask = scenarios$mask[iscen]

# iscen = as.integer(Sys.getenv("SGE_TASK_ID"))
# if (is.na(iscen)) stop("no scenario denoted")


outdir <- file.path(basedir, "Coregistered")
if (!file.exists(outdir)){
  dir.create(outdir, showWarnings=FALSE)
}

addstub = ""
if (mask) {
  addstub <- paste0(addstub, "_skullmask")
}

fnames = df$img
if (mask){
  fnames = df$ss
}

ifname <- fnames[1]
ref.img <- fnames[1]

  for (iname in 2:length(fnames)){
    ifname <- fnames[iname]

    stub <- paste0("rigid_", basename(ifname))
    ## strip off .nii.gz or .nii
    
    outfile <- file.path(outdir, paste0(nii.stub(stub), addstub))
    outmat = paste0(nii.stub(outfile), ".mat")

    opts = "-v"
    if (mask){
      opts = paste(opts, "-cost leastsq")
    }
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
