rm(list=ls())
library(oro.dicom)
library(oro.nifti)
library(plyr)
library(scales)

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
  study = "Registration_ICES"


  rootdir <- path.expand(rootdir)
  homedir <- file.path(rootdir, study)
  homedir <- path.expand(homedir)

#progdir <- file.path(dirname(basedir), "programs")
  progdir <- file.path(rootdir, "programs")
  source(file.path(progdir, "convert_DICOM.R"))
  source(file.path(progdir, "fslhd.R"))

setwd(homedir)
niis = list.files(path=homedir, 
  full.names=TRUE, 
  pattern="\\.nii\\.gz$",
  recursive=TRUE)

niis = niis[!grepl("Skull_Stripped|dcm2nii", niis)]

tars = list.files(path=homedir, 
  full.names=FALSE, 
  pattern="\\.tar\\.gz$",
  recursive=TRUE)

df = data.frame(tar=tars, stringsAsFactors=FALSE)
df$id = gsub("(.*)/Sorted/.*", "\\1", df$tar)
df$fname = basename(df$tar)
df$reg = gsub("_ungantry", "", df$fname)
df$gantry = grepl("_ungantry", df$fname)

drop.gantry = df$reg[df$gantry]

df = df[ !(df$fname %in% drop.gantry), ]

first = function(x) {
  if (NCOL(x) > 1) {
    x[1,]
  } else {
    x[1]
  }
}

firsts = ddply(df, .(id), first)

script = laply(firsts$tar, function(x){
  start = paste0("rsync -rav $jenig:$dex/Registration/", x, " ./")
})

writeLines(script, file.path(progdir, "Copy_First.sh"))