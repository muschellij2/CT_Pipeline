
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
  fromstudy = "MISTIE"
  fromdir = file.path(rootdir, fromstudy)
  
  # dropids = "100-318"
  dropids = NULL
  # keepids = c("100-318", "102-360")
  keepids = NULL

  getfiles <- function(basedir){
    files <- dir(path=basedir, pattern="dcm$|DCM$", full.names=TRUE, 
          recursive=TRUE)
    paths <- unique(dirname(files))
    return(list(files=files, paths=paths))
  }

  dirs <- list.dirs(fromdir, recursive=FALSE, full.names=FALSE)
  ptpat <- "\\d\\d\\d-(\\d|)\\d\\d\\d"
  ids <- grep(paste0(ptpat, "$"), dirs, value=TRUE)
  if (!is.null(dropids)) {
    ids = ids[ !(ids %in% dropids) ]
  }
  if (!is.null(keepids)) {
    ids = ids[ (ids %in% keepids) ]
  }  
  idir = 1;
  # for (idir in seq_along(ids)){
    files = getfiles(fromdir)

    xfiles = files$files
    if (!is.null(dropids)) {
      xfiles = xfiles[ !grepl(dropids, xfiles, fixed=TRUE) ]
    }
    if (!is.null(keepids)) {
      mat = sapply(keepids, grepl, x=xfiles)
      mat = matrix(mat, nrow=length(xfiles))
      mat = apply(mat, 1, any)
      xfiles = xfiles[ mat ]
    }    
    tofiles = gsub(fromstudy, "Registration", xfiles)
    paths = dirname(tofiles)
    make.dir = function(path){
      system(sprintf('mkdir -p "%s"', path))
    }
    # x = l_ply(unique(paths), make.dir, .progress="text")

    df = data.frame(xfiles, tofiles, stringsAsFactors=FALSE)
    cuts = as.numeric(cut(seq(nrow(df)), 5))
    df = df[ cuts == 4, ]
    # df = df[nrow(df):1,]
    x = mlply(df,function(xfiles,tofiles){
      file.copy(xfiles, tofiles)
    }, .progress = "text")
