rm(list=ls())
library(oro.dicom)
library(bitops)
library(arules)

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
  source(file.path(progdir, "convert_DICOM.R"))

  basedir <<- file.path(homedir, id)

}

#### setting up if things are on the cluster or not
# id <- "238-4136"
ids <- c("205-509", "205-517", "205-519", "225-502", "225-503", "225-504", 
"225-505", "225-506", "225-507", "225-510", "225-511", "225-515", 
"225-522", "225-523", "225-524", "232-512", "232-513", "232-514", 
"232-516", "232-521", "232-526", "289-518", "289-525", "301-520")

verbose=TRUE
untar = FALSE
convert <- TRUE
skullstrip <- TRUE
regantry <- TRUE
untgantry <- TRUE
runall <- FALSE

### initial setup
iid <- 2
id <- ids[iid]
setup(id)

zeroed <- dir(path=homedir, pattern= ".*Zeroed.*\\.nii\\.gz", recursive=TRUE, full.names=TRUE)
for (ifile in seq_along(zeroed)) system(sprintf('rm "%s"', zeroed[ifile]))

# zeroed <- dir(path=homedir, pattern= ":.*.gz", recursive=TRUE, full.names=TRUE)
# for (ifile in seq_along(zeroed)) system(sprintf('rm "%s"', zeroed[ifile]))

# zeroed <- dir(path=homedir, pattern= ":.*.txt", recursive=TRUE, full.names=TRUE)
# for (ifile in seq_along(zeroed)) system(sprintf('rm "%s"', zeroed[ifile]))

# zeroed <- dir(path=homedir, pattern= "'.*.tar.gz", recursive=TRUE, full.names=TRUE)
# for (ifile in seq_along(zeroed)) system(sprintf('rm "%s"', zeroed[ifile]))


if (regantry){
  ### re-run gantry tilt on the data
    files <- dir(path=homedir, pattern="_ungantry.tar.gz$", full.names=TRUE, 
        recursive=TRUE)
    if (untgantry) {
      ifile <- files[1]
      ### untarball the files using R's command
      if (length(files) > 0){
        for (ifile in files){
          dname <- dirname(ifile)
          untar(ifile, exdir = dname, compressed='gzip')
          if (verbose) print(ifile)
        }
      }
    } # untarball gantry
    fnames <- basename(files)  
    gantniis <- gsub("_ungantry.tar.gz", ".nii.gz", fnames, fixed=TRUE)
    gantniis <- file.path(dirname(dirname(files)), gantniis)
}

#### loop through IDS and convert them to nii, gantry tilted
### 301-520 needs to use Study Date/Time instead of Series Date/Time
for (iid in 1:length(ids)){
  id <- ids[iid]
  setup(id)
  # source(file.path(progdir, "file_functions.R"))
  dcmsortopt <- ifelse(id =="301-520", '-s ', "")


  infofile <- file.path(basedir, "Dropout_Information.Rda")
  file.remove(infofile)



  ### started 11:55
  contime <- NULL
  if (convert) {
    ### convert the dicoms
    contime <- system.time(convert_DICOM(basedir, progdir, 
                            verbose=verbose, untar=untar, dcmsortopt=dcmsortopt))

    ## dropout the niis that are not needed
    lis <- includeMatrix(basedir, dropstring="ungantry", error=TRUE)
    outs <- lis$outs
    mis <- lis$mis

    dropniis <- outs$fname[outs$Takeout]
    dropniis <- getBase(basename(dropniis), 1)

    if (length(dropniis) > 0){
      dropniis <- file.path(basedir, paste0(dropniis, ".nii.gz"))
      for (ifile in dropniis) system(sprintf('rm "%s"', ifile))
    }

    save(outs, mis, file = infofile)
  }
  
}


if (skullstrip){

  if (runall) {
    system.time(Skull_Strip(basedir, progdir, CTonly=TRUE, opts="-f 0.1", 
      verbose=verbose))
  } else if (regantry){
  ### skull strip the gantry tilt files only
    nniis <- length(gantniis)
    # nniis <- 300
    start <- 301
    for (ifile in start:nniis){
      nii <- gantniis[ifile]
      iddir <- dirname(dirname(nii))
      outdir <- file.path(iddir, "Skull_Stripped")
      Skull_Strip_file(img=nii, progdir=progdir, outdir=outdir, opts="-f 0.1", verbose=verbose)
    }
  } else {
    stop("Don't konw who to run - not runall and not just run gantry")
  }

}


# outs <- sapply(txts, getInfo)

# outs <- t(outs)
# rownames(outs) <- NULL
# outs <- data.frame(outs, stringsAsFactors=FALSE)
# outs$fname <- txts
# outs$SeriesNum <- as.numeric(outs$SeriesNum)
# outs$GantryDetectorTilt <- as.numeric(outs$GantryDetectorTilt)

# ## take out MRIs, localizers, dose reports, bone, cervical
# outs$Takeout <- FALSE
# outs$Takeout <- outs$Takeout | outs$Modality %in% "MR"
# outs$Takeout <- outs$Takeout | grepl("LOCALIZER", outs$itype)
# outs$Takeout <- outs$Takeout | grepl("SECONDARY", outs$itype)

# outs$Takeout <- outs$Takeout | grepl("Dose Report", outs$SeriesDesc)
# outs$Takeout <- outs$Takeout | grepl("BONE", outs$SeriesDesc)
# outs$Takeout <- outs$Takeout | grepl("Bone", outs$SeriesDesc)
# outs$Takeout <- outs$Takeout | grepl("CERVICAL", outs$StudyDesc)
# outs$Takeout <- outs$Takeout | grepl("Patient Protocol", outs$SeriesDesc)
# outs$Takeout <- outs$Takeout | grepl("SCOUT", outs$SeriesDesc)
# outs$Takeout <- outs$Takeout | grepl("CENTERING", outs$SeriesDesc)

# outs$Takeout <- outs$Takeout | outs$SeriesNum > 20

# outs$Takeout[is.na(outs$Takeout)] <- FALSE
# kept <- outs[! outs$Takeout, ]
# rownames(kept) <- NULL
# kept[ , c("itype", "SeriesDesc", "fname")]



# refdir="/Volumes/DATA/New_Age_Test/265-389"
# 
# cd $refdir
# for f in *.nii.gz; 
# do
# 


# sh ~/DHanley/CT_Registration/Brain_Seg_Rerun/programs/Brain_Seg_Function.sh -i "$f" -o ./Skull_Stripped
# 
# done;


# for (i in mis){
#   setwd(file.path(basedir, "Sorted"))
#   system(sprintf('tar -xf "%s.tar.gz"', i))
# }