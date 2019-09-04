rm(list=ls())
library(cttools)
library(fslr)
library(plyr)
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
ROIformat = FALSE
study = "MISTIE_III"
if (ROIformat) {
  study = paste0(study, "_ROI")
}

setup(study, study=study)

ids = list.dirs(homedir, recursive=FALSE, full.names=FALSE)
ids = basename(ids)
ids = grep("\\d\\d\\d-(\\d|)\\d\\d\\d", ids, value=TRUE)
length(ids)

verbose =TRUE
untar = FALSE
convert <- TRUE
skullstrip <- TRUE
plotss = FALSE
regantry <- FALSE
untgantry <- FALSE
runall <- TRUE
useRdcmsort= TRUE
useRdcm2nii= FALSE
removeDups = TRUE
isSorted = NULL
if (ROIformat) isSorted = FALSE
dcm2niicmd = "dcm2niix"

### initial setup
# iid <- length(ids)

iid <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if (is.na(iid)) iid <- 13

id <- ids[iid]
setup(id, study = study)

#### loop through IDS and convert them to nii, gantry tilted
### 301-520 needs to use Study Date/Time instead of Series Date/Time
# for (iid in 1:length(ids)){

setup(id, study=study)

# for (iid in seq_along(ids)){

  id <- ids[iid]
  print(id)
  setup(id, study = study)
  # source(file.path(progdir, "file_functions.R"))

  dcmsortopt <- ifelse(grepl("-301$", id), '-s ', "")
  useStudyDate <- grepl("-301$", id)


  if (convert) {
    ### convert the dicoms
  infofile <- file.path(basedir, "Dropout_Information.Rda")
  if (file.exists(infofile)) file.remove(infofile)
  
### time for convertsion
  contime <- NULL
  gf = getfiles(basedir)
# t = iconv("UTF-8","UTF-8//IGNORE",$t);
  
    if (length(gf$files) > 0 | untar){


      contime <- system.time(convert_DICOM(basedir, 
                              verbose=verbose, untar=untar, 
                              useRdcmsort= useRdcmsort, 
                              useRdcm2nii= useRdcm2nii,
                              id = id, 
                              isSorted = isSorted,
                              removeDups=removeDups,
                              dcmsortopt=dcmsortopt, 
                              ROIformat = ROIformat,
                              dcm2niicmd=dcm2niicmd,
                              gt_correct = FALSE, 
                              useStudyDate = useStudyDate,
                              change_transfer_syntax = TRUE))

 
      print(contime)


    # }
    }
  } ## end of if convert

# }


#### skull stripping
if (!ROIformat){
    # id <- ids[iid]
    # setup(id)
    if (skullstrip){

      if (runall) {
        cat("Skull Stripping")


        imgs = list.files(basedir, pattern = "\\.nii", 
          recursive=FALSE,
            full.names=TRUE)
        iimg = 1
        # int = c("0.01", "0.1", "0.35")
        int = "0.01"

        mid.folder = function(x, folname = ""){
          d = dirname(x)
          b = basename(x)
          file.path(d, folname, b)
        }

        presmooth=c(TRUE)
        refill = c(FALSE)

        for (iimg in seq_along(imgs)){
          
          img = imgs[iimg]

          ofile = nii.stub(img)
          ssdir = file.path(dirname(img), "Skull_Stripped")
          if (!file.exists(ssdir)){
            dir.create(ssdir)
          }
          ofile = mid.folder(ofile, "Skull_Stripped")
          ofile = paste0(ofile, "_SS")
                      
          app = "_nopresmooth"
          if (presmooth) app = ""


          re_app = ""
          if (refill) re_app = "_refill"

            
            outfile = paste0(ofile, "_", int, app, re_app)
            ext = get.imgext()
            if (!file.exists(paste0(outfile, ext))){
              x = CT_Skull_Strip(img = img, 
                                 outfile = outfile, 
                                 retimg=FALSE, verbose=TRUE, 
                                  # -w 2
                                 opts=paste0("-f ", int, " -v"),
                                 inskull_mesh = FALSE,
                                 refill = refill,
                                 presmooth=presmooth)
            } # file exists

        } # end iimg

      } else {

        outdir <- file.path(basedir, "Skull_Stripped")
        niis <- dir(path=basedir, pattern=".nii.gz", 
          full.names=TRUE, recursive=FALSE)

      }
    }


    ### make plots of overlays for skull stripping 
    if (plotss){

      ssdir = file.path(basedir, "Skull_Stripped")
      odir = file.path(ssdir, "overlays")
      dir.create(odir, showWarnings=FALSE)
      ss.imgs = list.files(ssdir, pattern=".*\\.nii\\.gz", 
        full.names=TRUE)
      imgs = gsub("_SS_(0\\.1|0\\.35|0\\.01)_Mask", "", ss.imgs, 
        fixed=FALSE)
      imgs = gsub("Skull_Stripped/", "", imgs)

      df = data.frame(img=imgs, img.mask =ss.imgs, 
        stringsAsFactors=FALSE)

      keep = apply(sapply(df, file.exists), 1, all)

      df = df[keep, ]

      # df = df[ grepl("0.01", df$img.mask, fixed=TRUE),]
      # img = imgs[1]
      # img.mask = ss.imgs[1]
      no.gz = function(x){
        x = gsub("\\.gz$", "", x)
        x = gsub("\\.nii$", "", x)
        x
      }
      
      m_ply(function(img, img.mask){
        fname = no.gz(basename(img.mask))
        ifname = no.gz(basename(img))

        outfile = file.path(odir, paste0(fname, ".png"))
        ff = gsub("_SS_Mask", "_SS", fname)
        ssfile = file.path(odir, paste0(ff, ".png"))
        imgfile = file.path(odir, paste0(ifname, ".png"))

        iimg = readNIfTI(img, reorient=FALSE)        
        iimg.mask = readNIfTI(img.mask, reorient=FALSE)        
        ssimg = iimg 
        ssimg[!(iimg.mask > 0)] = 0

        try({
          png(outfile, type="cairo")
           mask.overlay(iimg, iimg.mask)
          dev.off()

            iimg@cal_min = 0
            iimg@cal_max = 100
          png(imgfile, type="cairo")
           ortho2(iimg)
          dev.off()          
            # img[ img >= window[2] | img < window[1] ] = 0
            ssimg@cal_min = 0
            ssimg@cal_max = 100

          png(ssfile, type="cairo")
           ortho2(ssimg)
          dev.off()          
        })
      }, .data=df, .progress = "text")
      
    }

    # mask.overlay(imgs[1], ss.imgs[1])
 
}
