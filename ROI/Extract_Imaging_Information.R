##################
#### Get statistics on first 111 images
#################

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

df = ddply(df, .(id), function(x) x[1,])
stopifnot(nrow(df) == 111)
df = df[, c("pid", "iddir", "raw", "id")]
df$raw = gsub("\\.nii$", "_Header_Info.Rda", basename(df$raw))
df$rda = file.path(df$iddir, "Sorted", df$raw)

get.val = function(rda, val){
  if ("dcmtables" %in% ls()) rm(list="dcmtables")
  ungant.rda = gsub("_Header_Info\\.Rda", 
    "_ungantry_Header_Info\\.Rda", 
    rda)
  if (file.exists(ungant.rda)) rda = ungant.rda
  load(rda)
  cn = colnames(dcmtables)
  n.slice = length(unique(dcmtables[, "0018-0050-SliceThickness"]))
  co.kern = unique(dcmtables[, val])
  co.kern$n.slice = n.slice
  stopifnot(nrow(co.kern) == 1)
  co.kern
}
rda = df$rda[1]
data = ldply(.data=df$rda,  get.val, 
  val=c("0018-1210-ConvolutionKernel", 
  "0008-0070-Manufacturer", 
  "0018-1120-GantryDetectorTilt"), 
.progress="text")


colnames(data) = c("kern", "man", "tilt", "nslices")
data$tilt = as.numeric(data$tilt)
# data$rda = df$rda
data$id = as.numeric(gsub("-", "", df$id))

write.csv(data, file.path(rootdir, "Registration", 
  "Imaging_Information.csv"), row.names=FALSE)
