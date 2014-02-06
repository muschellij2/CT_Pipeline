rm(list=ls())
library(oro.dicom)
library(oro.nifti)
library(plyr)
library(scales)
library(reshape2)

#### delete all ROI files
### find . -regextype posix-extended -regex "^./[0-9].*[0-9]$"
###  -exec rm -r {} \;

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

#progdir <- file.path(dirname(basedir), "programs")
  progdir <- file.path(rootdir, "programs")
  source(file.path(progdir, "convert_DICOM.R"))
  source(file.path(progdir, "fslhd.R"))


#### setting up if things are on the cluster or not

load(file.path(rootdir, "Registration", 
  "Registration_Image_Names.Rda"))
xdf = df
reorient = FALSE
normalize = FALSE
subsamp = FALSE
imgs = mlply(.fun = function(outfile, roi.nii, raw, ss){
  c(outfile, roi.nii, raw, ss)
}, .data=df[, c("outfile", "roi.nii", "raw", "ss")])
attr(imgs, "split_labels") <- NULL

x = imgs[[2]]
if (reorient){
  rets = laply(.data=imgs, .fun = function(x){
    acpc_reorient(infiles = x)
  }, .progress="text")
  stopifnot(all(rets == 0))
  rois = df$roi.nii
  l_ply(.data=rois, .fun = fslthresh, unzip=TRUE, 
    .progress = "text")
  # acpc_reorient(infiles = x)
}


imgs = mlply(.fun = function( raw, roi.nii){
  c(raw, roi.nii)
}, .data=df[, c("raw", "roi.nii")])
attr(imgs, "split_labels") <- NULL

x = imgs[[2]]

if (normalize){
  ### matlab problems if I don't change directories
  gwd = getwd()
  setwd("~/")
  rets = laply(.data=imgs, .fun = function(x){
      r = run_ctnorm(rawfile = x[1], roifile = x[2])
      return(r)
    }, .progress="text")
  setwd(gwd)
}

### after orientation

##### subsampling to 2mm - for atlas overlap
df = xdf[, c("roi.nii", "raw")]
df$roi.nii = file.path(dirname(df$roi.nii), 
  paste0("bws", basename(df$roi.nii)))
df$raw = file.path(dirname(df$raw), 
  paste0("w", basename(df$raw)))

exists = t(apply(df, 1, file.exists))
aexists = apply(exists, 1, all)
nonorm = xdf[which(!aexists),]
nonorm = gsub(rootdir, "", nonorm$copydir)
nonorm = gsub("reoriented", "", nonorm)

sum(!aexists)
df = df[aexists,]
df$id = 1:nrow(df)
melted = melt(df, id.vars="id")

melted$outfile = file.path(dirname(melted$value), 
  paste0("2mm_", basename(melted$value)))

if (subsamp){

    value = melted$value[1]
    outfile = melted$outfile[1]
    m_ply(melted[, c("value", "outfile")], 
      .fun = function(value, outfile){
        # outdir = dirname(outfile)
        # niigz = list.files(outdir, full.names=TRUE, 
        # pattern="*.nii.gz")
        # file.remove(niigz)
        file.remove(outfile)
        fslsub2(value, outfile, unzip=TRUE)
      }, 
    .progress= "text")
}


