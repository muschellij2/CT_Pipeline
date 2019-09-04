rm(list=ls())
library(cttools)
library(fslr)
library(plyr)
library(methods)
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
## ROIformat after 134-327.zip
ROIformat = TRUE
study = "Registration"
if (ROIformat) {
  study = "ROI_data"
}

setup(study, study=study)

if (ROIformat) setup("Long_ROIS", study=study)


ids = list.dirs(homedir, recursive=FALSE, full.names=FALSE)
ids = basename(ids)
ids = grep("\\d\\d\\d-(\\d|)\\d\\d\\d", ids, value=TRUE)
length(ids)



iid <- as.numeric(Sys.getenv("SGE_TASK_ID"))

# iid = 100 needs ss

if (is.na(iid)) iid <- 118


all.df = NULL
# for (iid in seq_along(ids)){
for (iid in seq_along(ids)){
  
  id <- ids[iid]
  print(id)
  setup(id, study = study)
  # source(file.path(progdir, "file_functions.R"))


  iddir = file.path(rootdir, "Registration", id)
  roidir = file.path(rootdir, "ROI_data", id)

  id.outdir = file.path(roidir, "overlays")

  if (!file.exists(id.outdir)){
    dir.create(id.outdir)
  }

  rois = list.files(roidir, pattern=".nii.gz", recursive=FALSE)
  imgs = file.path(iddir, gsub("ROI.nii", ".nii", rois))
  rois = file.path(roidir, rois)

  if (length(rois) > 0){
    df = data.frame(img=imgs, roi = rois, 
      iddir = iddir, roidir = roidir, stringsAsFactors=FALSE)

    bad.df = df[!file.exists(df$img) & file.exists(df$roi),]
    if (nrow(bad.df) > 0){
      print(bad.df)
    }


    df = df[ file.exists(df$img) & file.exists(df$roi), ]

    irow = 1

    if (nrow(df) > 0) {
        df$good = FALSE
      for (irow in seq(nrow(df))){

        pngname = file.path(id.outdir, 
          paste0(nii.stub(df$roi[irow], bn=TRUE),
            ".png"))

        roi.f = df$roi[irow]
        img.f = df$img[irow]

        roi.dims = as.numeric(sapply(c("dim1", "dim2", "dim3"), 
          function(x) {
            fslval(roi.f, x, verbose=FALSE)
        }))
        img.dims = as.numeric(sapply(c("dim1", "dim2", "dim3"), 
          function(x) {
            fslval(img.f, x, verbose=FALSE)
        }))      

        dim.ok = all(roi.dims == img.dims)

        if ( (!file.exists(pngname)) & dim.ok){
          roi = readNIfTI(roi.f, reorient=FALSE)
          img = readNIfTI(img.f, reorient=FALSE) 
          xyz <- ceiling(dim(img)/2)
          if (sum(roi) > 0){
            xyz=cog(roi, ceil=TRUE)
          }
          png(pngname, type="cairo")
            mask.overlay(img, roi, window=c(0, 100), xyz=xyz)
          dev.off()
          print(roi.f)
        } else {
          if (!file.exists(pngname)){
            print(paste0("Bad region ", roi.f))
            print(roi.dims)
            print(img.dims)
          }
        }

        if (file.exists(pngname)){
          df$good[irow] = TRUE
        }

      } # end irow
      all.df = rbind(all.df, df)
    } # end if df > 0
  } # end if length rois > 0
}

all.df$png = file.path(all.df$roidir, "overlays", 
  paste0(nii.stub(all.df$roi, bn=TRUE), ".png"))

all.df$id = basename(all.df$iddir)
resdir = file.path(rootdir, "Registration", "results")
save(all.df, 
  file=file.path(resdir, "ROI_Filenames.Rda")
  )


