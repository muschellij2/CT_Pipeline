rm(list=ls())
library(oro.dicom)
library(oro.nifti)
library(plyr)
library(scales)

study = "Registration_ICES"


setup <- function(id, study = "Registration"){
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
  homedir <<- file.path(rootdir, study)
  homedir <<- path.expand(homedir)

#progdir <- file.path(dirname(basedir), "programs")
  progdir <<- file.path(rootdir, "programs")
  source(file.path(progdir, "convert_DICOM.R"))
  source(file.path(progdir, "fslhd.R"))

  basedir <<- file.path(homedir, id)

}


#### setting up if things are on the cluster or not
# id <- "238-4136"
# ids <- c("205-509", "205-517", "205-519", "225-502", "225-503", "225-504", 
# "225-505", "225-506", "225-507", "225-510", "225-511", "225-515", 
# "225-522", "225-523", "225-524", "232-512", "232-513", "232-514", 
# "232-516", "232-521", "232-526", "289-518", "289-525", "301-520")
if (study %in% "Registration") {
  ids = all.ids = c("100-318", "100-362", "100-365", "101-306", "101-307", "101-308", 
"102-317", "102-322", "102-323", "102-324", "102-326", "102-331", 
"102-347", "102-349", "102-351", "102-360", "102-367", "102-374", 
"102-391", "102-393", "102-403", "102-406", "111-415", "120-376", 
"131-310", "131-316", "131-334", "131-354", "133-409", "133-417", 
"134-304", "134-305", "134-320", "134-327", "134-338", "134-343", 
"134-345", "134-379", "134-380", "134-381", "134-382", "134-392", 
"134-408", "134-411", "134-412", "134-416", "152-302", "152-303", 
"152-348", "152-353", "157-328", "157-329", "157-332", "157-335", 
"157-336", "157-370", "157-372", "157-399", "157-410", "161-413", 
"173-312", "173-313", "173-325", "173-339", "173-341", "173-361", 
"173-364", "173-368", "173-384", "173-396", "173-404", "175-387", 
"175-397", "175-405", "179-359", "179-373", "179-383", "179-386", 
"179-394", "179-395", "179-401", "179-402", "184-342", "184-388", 
"191-301", "191-309", "191-311", "191-314", "191-315", "191-319", 
"191-321", "191-330", "191-333", "191-375", "191-378", "191-400", 
"210-344", "216-390", "216-414", "219-340", "219-350", "222-337", 
"222-357", "222-358", "223-355", "223-369", "223-407", "230-346", 
"230-352", "230-356", "230-363", "230-366", "230-371", "230-377", 
"234-385", "265-389", "265-398")
} else {
  ids = all.ids = c("205-509", "205-517", "205-519", "225-502", "225-503", "225-504", 
"225-505", "225-506", "225-507", "225-510", "225-511", "225-515", 
"225-522", "225-523", "225-524", "232-512", "232-513", "232-514", 
"232-516", "232-521", "232-526", "289-518", "289-525", "301-520")
}

# bad.ids = c(16, 27, 28, 33, 38, 40, 43, 57, 59, 61, 62, 63, 64, 65, 66, 
# 67, 68, 69, 71, 72, 75, 76, 77, 78, 79, 80, 81, 82, 85, 98, 99, 
# 100, 105, 107, 108, 109, 110, 111, 112, 113, 114, 115)
bad.ids = c(21, 47, 48, 49, 5, 50, 54, 60, 69, 84, 97, 98)
# ids = "100-362"


runonlybad = FALSE

if (runonlybad) ids = ids[bad.ids]

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
dcm2niicmd = "dcm2nii_2009"

### initial setup
# iid <- length(ids)

iid <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if (is.na(iid)) iid <- 19

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
          untar(ifile, exdir = dname, compressed='gzip')
          if (verbose) print(ifile)
        }
      }
    } # untarball gantry
    fnames <- basename(files)
    gantniis <- gsub("_ungantry.tar.gz", ".nii.gz", 
      fnames, fixed=TRUE)
    gantniis <- file.path(dirname(dirname(files)), gantniis)
} 

#### loop through IDS and convert them to nii, gantry tilted
### 301-520 needs to use Study Date/Time instead of Series Date/Time
# for (iid in 1:length(ids)){

  id <- ids[iid]
  setup(id, study = study)
  # source(file.path(progdir, "file_functions.R"))
  dcmsortopt <- ifelse(id =="301-520", '-s ', "")




  if (convert) {
    ### convert the dicoms
  infofile <- file.path(basedir, "Dropout_Information.Rda")
  file.remove(infofile)
  
### time for convertsion
  contime <- NULL

    contime <- system.time(convert_DICOM(basedir, progdir, 
                            verbose=verbose, untar=untar, 
                            useRdcmsort= useRdcmsort, 
                            useRdcm2nii= useRdcm2nii,
                            id = id, 
                            removeDups=removeDups,
                            dcmsortopt=dcmsortopt, 
                            dcm2niicmd=dcm2niicmd))

    # contime <- system.time(convert_DICOM(basedir, progdir, 
    #                         verbose=verbose, untar=untar, 
    #                         useRdcmsort= TRUE, 
    #                         useRdcm2nii= FALSE,
    #                         id = id, 
    #                         dcmsortopt=dcmsortopt, 
    #                         dcm2niicmd=dcm2niicmd))    
    print(contime)

    ## dropout the niis that are not needed
    # lis <- includeMatrix(basedir, dropstring="ungantry", error=TRUE)
    # outs <- lis$outs
    # mis <- lis$mis

    # dropniis <- outs$fname[outs$Takeout]
    # dropniis <- getBase(basename(dropniis), 1)

    # if (length(dropniis) > 0){
    #   dropniis <- file.path(basedir, paste0(dropniis, ".nii.gz"))
    #   for (ifile in dropniis) system(sprintf('rm "%s"', ifile))
    # }

    # save(outs, mis, file = infofile)
  # }
  
}


#### skull stripping
# for (iid in 1:length(ids)){
  # id <- ids[iid]
  # setup(id)
  if (skullstrip){

    if (runall) {

      system.time(Skull_Strip(basedir, progdir, CTonly=TRUE, 
        opts="-f 0.1 -b", 
        verbose=verbose))


      system.time(Skull_Strip(basedir, progdir, CTonly=TRUE, 
        opts="-f 0.01 -b", 
        verbose=verbose))

      system.time(Skull_Strip(basedir, progdir, CTonly=TRUE, 
        opts="-f 0.35 -b", 
        verbose=verbose))


    }
  }
# }


### make plots of overlays for skull stripping 
if (plotss){

  ssdir = file.path(basedir, "Skull_Stripped")
  odir = file.path(ssdir, "overlays")
  dir.create(odir, showWarnings=FALSE)
  ss.imgs = list.files(ssdir, pattern=".*\\.nii\\.gz", 
    full.names=TRUE)
  imgs = gsub("_SS_Mask_(0\\.1|0\\.35|0\\.01)", "", ss.imgs, fixed=FALSE)
  imgs = gsub("Skull_Stripped/", "", imgs)

  df = data.frame(img=imgs, img.mask =ss.imgs, 
    stringsAsFactors=FALSE)

  keep = apply(sapply(df, file.exists), 1, all)

  df = df[keep, ]

  # df = df[ grepl("0.01", df$img.mask, fixed=TRUE),]
  # img = imgs[1]
  # img.mask = ss.imgs[1]
  m_ply(function(img, img.mask){
    fname = gsub("\\.gz", "", basename(img.mask))
    fname = gsub("\\.nii$", "", fname)
    outfile = file.path(odir, paste0(fname, ".png"))
    png(outfile, type="cairo")
      mask.overlay(img, img.mask)
    dev.off()
  }, .data=df, .progress = "text")
  
}

  # mask.overlay(imgs[1], ss.imgs[1])
 


# if (skullstrip){

#   if (regantry){
#   ### skull strip the gantry tilt files only
#     nniis <- length(gantniis)
#     # nniis <- 300
#     start <- 1
#     for (ifile in start:nniis){
#       nii <- gantniis[ifile]
#       iddir <- dirname(dirname(nii))
#       outdir <- file.path(iddir, "Skull_Stripped")
#       Skull_Strip_file(img=nii, progdir=progdir, outdir=outdir, 
# opts="-f 0.1", verbose=verbose)
#     }
#   } else {
#     stop("Don't konw who to run - not runall and not just run gantry")
#   }

# }


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


# sh ~/DHanley/CT_Registration/Brain_Seg_Rerun/programs/Brain_Seg_Function.sh \
#  -i "$f" -o ./Skull_Stripped
# 
# done;

# SSPLOT.e1616337.3: in if (!all(dim(x)[1:3] == dim(y)[1:3])) { : 
# SSPLOT.e1616337.68: in if (!all(dim(x)[1:3] == dim(y)[1:3])) { : 
# SSPLOT.e1616337.75: in if (!all(dim(x)[1:3] == dim(y)[1:3])) { : 


# for (i in mis){
#   setwd(file.path(basedir, "Sorted"))
#   system(sprintf('tar -xf "%s.tar.gz"', i))
# }