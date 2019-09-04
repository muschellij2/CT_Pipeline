rm(list=ls())
library(cttools)
library(fslr)
library(plyr)
library(methods)
options(matlab.path='/Applications/MATLAB_R2016a.app/bin')


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
ROIformat = FALSE
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


# runonlybad = FALSE

# if (runonlybad) ids = ids[bad.ids]

verbose =TRUE
untar = FALSE
convert <- TRUE
skullstrip <- TRUE
plotss = TRUE
regantry <- FALSE
untgantry <- FALSE
runall <- TRUE
useRdcmsort= TRUE
useRdcm2nii= FALSE
removeDups = TRUE
isSorted = NULL
if (ROIformat) isSorted = FALSE
dcm2niicmd = "dcm2nii_2009"

### initial setup
# iid <- length(ids)

iid <- as.numeric(Sys.getenv("SGE_TASK_ID"))

# iid = 100 needs ss

if (is.na(iid)) iid <- 51
id <- ids[iid]
setup(id, study = study)

zeroed <- dir(path=basedir, pattern= ".*Zeroed.*\\.nii\\.gz", 
  recursive=TRUE, full.names=TRUE)
for (ifile in seq_along(zeroed)) {
  system(sprintf('rm "%s"', zeroed[ifile]))
}

# zeroed <- dir(path=homedir, pattern= ":.*.gz", 
# recursive=TRUE, full.names=TRUE)
# for (ifile in seq_along(zeroed)) system(sprintf('rm "%s"', 
# zeroed[ifile]))

# zeroed <- dir(path=homedir, pattern= ":.*.txt", recursive=TRUE, 
# full.names=TRUE)
# for (ifile in seq_along(zeroed)) system(sprintf('rm "%s"', 
# zeroed[ifile]))

# zeroed <- dir(path=homedir, pattern= "'.*.tar.gz", 
# recursive=TRUE, full.names=TRUE)
# for (ifile in seq_along(zeroed)) system(sprintf('rm "%s"', 
# zeroed[ifile]))


if (regantry){
  ### re-run gantry tilt on the data
   files <- dir(path=homedir, pattern="_ungantry.tar.gz$", 
    full.names=TRUE, 
        recursive=TRUE)
    if (untgantry) {
      ifile <- files[1]
      ### untarball the files using R's command
      if (length(files) > 0){
        for (ifile in files){
          dname <- dirname(ifile)
          untar(ifile, exdir = dname, 
            compressed="gzip")
          if (verbose) print(ifile)
        }
      }
    } # untarball gantry
    fnames <- basename(files)
    gantniis <- gsub("_ungantry.tar.gz", ".nii.gz", 
      fnames, fixed=TRUE)
    gantniis <- file.path(dirname(dirname(files)), gantniis)
} 

## loop through IDS and convert them to nii, gantry tilted
## 301-520 needs to use Study Date/Time instead of Series Date/Time
# for (iid in 1:length(ids)){

setup(id, study=study)

# for (iid in 60:70){
  id <- ids[iid]
  print(id)
  setup(id, study = study)
  # source(file.path(progdir, "file_functions.R"))
  dcmsortopt <- ifelse(id %in% c("301-520", "191-309"), 
    '-s ', "")
  useStudyDate <- id %in% c("301-520", "191-309")




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
                              removeDups = removeDups,
                              dcmsortopt = dcmsortopt, 
                              ROIformat = ROIformat,
                              dcm2niicmd = dcm2niicmd, 
                              useStudyDate = useStudyDate, 
                              check_series = FALSE,
                              useNA='ifany', 
                              add.img.dir = TRUE))


      # contime <- system.time(convert_DICOM(basedir, progdir, 
      #                         verbose=verbose, untar=untar, 
      #                         useRdcmsort= TRUE, 
      #                         useRdcm2nii= FALSE,
      #                         id = id, 
      #                         dcmsortopt=dcmsortopt, 
      #                         dcm2niicmd=dcm2niicmd))    
      print(contime)

      ## dropout the niis that are not needed
      # lis <- includeMatrix(basedir, dropstring="ungantry", 
      # error=TRUE)
      # outs <- lis$outs
      # mis <- lis$mis

      # dropniis <- outs$fname[outs$Takeout]
      # dropniis <- getBase(basename(dropniis), 1)

      # if (length(dropniis) > 0){
      #  dropniis <- file.path(basedir, paste0(dropniis, ".nii.gz"))
      #  for (ifile in dropniis) system(sprintf('rm "%s"', ifile))
      # }

      # save(outs, mis, file = infofile)
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

        cat("Skull Stripping\n")
          # system.time(Skull_Strip(basedir, CTonly=TRUE, 
          #   opts="-f 0.1 -b", 
          #   verbose=verbose))


          # system.time(Skull_Strip(basedir, CTonly=TRUE, 
          #   opts="-f 0.01 -b", 
          #   verbose=verbose))

          # system.time(Skull_Strip(basedir, CTonly=TRUE, 
          #   opts="-f 0.35 -b", 
          #   verbose=verbose))

        imgs = list.files(basedir, pattern = "\\.nii", 
          recursive=FALSE,
          full.names=TRUE)
        iimg = 1
        
        scen = expand.grid(int=c("0.01", "0.1", "0.35"),
                           presmooth=c(TRUE, FALSE),
                           refill = c(FALSE))
        rownames(scen)= NULL
        w = !scen$presmooth & scen$refill
        scen = scen[!w, ]

        mid.folder = function(x, folname = ""){
          d = dirname(x)
          b = basename(x)
          file.path(d, folname, b)
        }

        ssdir = file.path(basedir, "Skull_Stripped")
        if (!file.exists(ssdir)){
          dir.create(ssdir)
        }
        for (iimg in seq_along(imgs)){
          
          img = imgs[iimg]

          ofile = nii.stub(img)
          ofile = mid.folder(ofile, "Skull_Stripped")
          ofile = paste0(ofile, "_SS")
          
          iscen = 1
          
          for (iscen in seq(nrow(scen))){
            
            int = scen$int[iscen]
            presmooth = scen$presmooth[iscen]
            refill = scen$refill[iscen]
            
            app = "_nopresmooth"
            if (presmooth) app = ""


            re_app = ""
            if (refill) re_app = "_refill"

            
            outfile = paste0(ofile, "_", int, app, re_app)
            x = CT_Skull_Strip(img = img, 
                               outfile = outfile, 
                               retimg=FALSE, verbose=TRUE, 
                               opts=paste0("-f ", int, " -v"),
                               inskull_mesh = FALSE,
                               refill = refill,
                               refill.thresh = .75,
                               presmooth=presmooth)   
            

          } # scen

        }
      } else {

        outdir <- file.path(basedir, "Skull_Stripped")
        niis <- dir(path=basedir, pattern=".nii.gz", 
          full.names=TRUE, recursive=FALSE)

        # # for (inii in niis){
        # inii = niis[1]
        #   Skull_Strip_file(img=inii,  
        #     outdir=outdir, opts="-f 0.1 -b", verbose=verbose)
        #   Skull_Strip_file(img=inii,
        #     outdir=outdir, opts="-f 0.01 -b", verbose=verbose)       
        #   Skull_Strip_file(img=inii, 
        #     outdir=outdir, opts="-f 0.35 -b", verbose=verbose)  
        # # }

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

# grep Err RERUN.o1800157.* | sed 's/.*157.\([0-9]\{1,3\}\):.*/\1/'
# grep Fail RERUN.o1800157.* | sed 's/.*157.\(.*\):\[.*/\1/'