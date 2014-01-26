rm(list=ls())
library(oro.dicom)
library(oro.nifti)
library(plyr)

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

  study = "Registration"
  todir <- file.path(rootdir, study)
  fromdir = file.path(rootdir, "MISTIE")
  
  dropids = "100-318"



  getfiles <- function(basedir){
    files <- dir(path=basedir, pattern="dcm$|DCM$", full.names=TRUE, 
          recursive=TRUE)
    paths <- unique(dirname(files))
    return(list(files=files, paths=paths))
  }

  dirs <- list.dirs(fromdir, recursive=FALSE, full.names=FALSE)
  ptpat <- "\\d\\d\\d-(\\d|)\\d\\d\\d"
  ids <- grep(paste0(ptpat, "$"), dirs, value=TRUE)
  ids = ids[ !(ids %in% dropids) ]

  idir = 1;
  # for (idir in seq_along(ids)){
    files = getfiles(fromdir)

    xfiles = files$files
    xfiles = xfiles[ !grepl(dropids, xfiles, fixed=TRUE) ]
    
    tofiles = gsub("MISTIE", "Registration", xfiles)
    paths = dirname(tofiles)
    make.dir = function(path){
      system(sprintf('mkdir -p "%s"', path))
    }
    # l_ply(unique(paths), make.dir, .progress="text")

    df = data.frame(xfiles, tofiles, stringsAsFactors=FALSE)
    x = mlply(df,function(xfiles,tofiles){
      file.copy(xfiles, tofiles)
    }, .progress = "text")
