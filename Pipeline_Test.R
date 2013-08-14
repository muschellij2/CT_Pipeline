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
id <- "301-520"
setup(id)
# source(file.path(progdir, "file_functions.R"))


infofile <- file.path(basedir, "Dropout_Information.Rda")
file.remove(infofile)

verbose=TRUE
untar = FALSE
convert <- TRUE
skullstrip <- TRUE
dcmsortopt <- '-s'

### started 11:55
contime <- NULL
if (convert) {
  contime <- system.time(convert_DICOM(basedir, progdir, 
                          verbose=verbose, untar=untar, dcmsortopt=dcmsortopt))


lis <- includeMatrix(basedir, dropstring="ungantry")
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

# ss <- strsplit(outs$itype, "\\\\")



print(contime)
# , dropstring = c("_CTA_")
if (skullstrip) system.time(Skull_Strip(basedir, progdir, CTonly=TRUE, opts="-f 0.1", 
  verbose=verbose))


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