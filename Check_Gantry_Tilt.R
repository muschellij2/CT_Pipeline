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
regantry <- FALSE
untgantry <- TRUE
runall <- FALSE

### initial setup
iid <- 1
id <- ids[iid]
setup(id)


all.outs <- NULL
#### loop through IDS and convert them to nii, gantry tilted
### 301-520 needs to use Study Date/Time instead of Series Date/Time
for (iid in 1:length(ids)){
  id <- ids[iid]
  setup(id)
  # source(file.path(progdir, "file_functions.R"))
  sortdir <- file.path(basedir, "Sorted")
  ## dropout the niis that are not needed

  infofile <- file.path(basedir, "Dropout_Information.Rda")
  load(infofile)

  outs <- outs[ outs$Modality %in% "CT", ]


  outs$file <- outs$fname
  outs$fname <- basename(outs$file)
  outs$dir <- dirname(outs$file)
  ss <- strsplit(outs$fname, split="_")

  gantry <- gsub("_ungantry", "", mis)

  gantry <- data.frame(gtilt=gantry, stringsAsFactors=FALSE)

  outs$DT <- sapply(ss, function(x) paste0(x[1:5], sep="", collapse="_"))

### get assessment of gantry tilt of true data
  tab <- table(DT=outs$DT)
  stopifnot(all(tab %in% c(1,2)))
  dt <- as.data.frame(tab)
  outs <- merge(outs, dt, by="DT", all.x=TRUE, sort=FALSE)
  
  drop <- (!grepl("_ungantry", outs$fname)) & outs$Freq == 2
  outs <- outs[ !drop, ]
  outs$Freq <- NULL

### take out recon data
  ss <- strsplit(outs$fname, split="_")
  outs$DT <- sapply(ss, function(x) paste0(x[1:4], sep="", collapse="_"))

  tab <- table(DT=outs$DT)
  dt <- as.data.frame(tab)
  outs <- merge(outs, dt, by="DT", all.x=TRUE, sort=FALSE)

  outs$SD <- tolower(outs$SeriesDesc)
  drop <- outs$Freq > 1 & grepl("recon", outs$SD)
  outs <- outs[ !drop, ] 
  outs$SD <- outs$Freq <- NULL

## take out any non-readable or whatever CT Data
  #outs <- outs[! outs$Takeout, ]
  rownames(outs) <- NULL

  ss <- strsplit(outs$fname, split="_")
  outs$ID <- sapply(ss, function(x) x[1])

  all.outs <- rbind(all.outs, outs)
  cat(id, "\n")
}


xall.outs <- all.outs


all.outs <- xall.outs[!xall.outs$Takeout, ] 

save(all.outs, file=file.path(homedir, "results", "Readable_Scans.rda"))


tilts <- all.outs$GantryDetectorTilt
tab <- table(tilts, useNA='ifany')
mean(tilts != 0)
sum(tilts != 0)
range(tilts)
nrow(all.outs)
png(file.path(resdir, "Gantry_Tilt_Distribution.png"))
hist(tilts[tilts!=0], main="Histogram of non-zero tilts", xlab="Gantry Tilt (Degrees)")
dev.off()

