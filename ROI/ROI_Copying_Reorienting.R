rm(list=ls())
library(fslr)
library(cttools)
library(plyr)
library(scales)
library(reshape2)
library(R.utils)
#### delete all ROI files
### find . -regextype posix-extended -regex "^./[0-9].*[0-9]$"
###  -exec rm -r {} \;

setup <- function(id, study="ROI_data"){
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

  basedir <<- file.path(rootdir, "Registration")

#progdir <- file.path(dirname(basedir), "programs")
  progdir <<- file.path(rootdir, "programs")


}

rerun = TRUE
#### setting up if things are on the cluster or not
verbose =TRUE
untar = TRUE
convert <- TRUE
skullstrip <- FALSE
plotss = TRUE
regantry <- FALSE
untgantry <- FALSE
runall <- TRUE
useRdcmsort= TRUE
useRdcm2nii= FALSE
removeDups = TRUE
ROIformat = TRUE
dcm2niicmd = "dcm2nii_2009"

### initial setup
# iid <- length(ids)

#### loop through IDS and convert them to nii, gantry tilted
### 301-520 needs to use Study Date/Time instead of Series Date/Time
# for (iid in 1:length(ids)){
study= "ROI_data"
setup("ROIS", study=study)

ids = list.dirs(homedir, recursive=FALSE, full.names=FALSE)
ids = grep("\\d\\d\\d-(\\d|)\\d\\d\\d", ids, value=TRUE)
length(ids)

# iid <- as.numeric(Sys.getenv("SGE_TASK_ID"))

# if (is.na(iid)) iid <- 1

# # for (iid in seq_along(ids)){

#   id = ids[iid]

#   setup(id)

niis = list.files(path = homedir, 
  pattern= "\\.nii\\.gz", 
  full.names=TRUE, recursive=TRUE)

niis = niis[!grepl("Skull_Stripped|dcm2nii", niis)]

bn = basename(niis)
dn = dirname(niis)

id.vals = basename(dn)

pid = as.numeric(gsub("(\\d\\d\\d)-((\\d|)\\d\\d\\d)", 
  "\\2", id.vals))
# topdirs = sapply(pid, function(x) ifelse(x < 500, "Registration", 
#   "Registration_ICES"))
topdirs = sapply(pid, function(x) "Registration")
raw.nii = gsub("ROI\\.nii", ".nii", bn)

iddir = file.path(rootdir, topdirs, id.vals)
copydir = file.path(iddir, "reoriented")

make.dir = function(path){
  system(sprintf('mkdir -p "%s"', path))
}
l_ply(unique(copydir), make.dir, .progress="text")
df = data.frame(roi.nii = niis, raw = raw.nii, copydir, pid, 
  id = id.vals, 
  iddir = iddir,
  stringsAsFactors=FALSE)

df$ss = gsub("\\.nii\\.gz", "_SS_0.01_Mask.nii.gz", df$raw)
df$ss = file.path(df$iddir, "Skull_Stripped", df$ss)
df$raw = file.path(df$iddir, df$raw)


exists = t(apply(df[, c("roi.nii", "raw", "ss")], 1, file.exists))
aexists = apply(exists, 1, all)
# stopifnot(all(aexists))
df[which(!aexists),]
df = df[aexists, ]


#### save names of .nii files
xdf = df

df$raw = file.path(df$copydir, basename(df$raw))
df$ss = file.path(df$copydir, basename(df$ss))
df$roi.nii = file.path(df$copydir, basename(df$roi.nii))


df$raw = gsub("\\.gz$", "", df$raw)
df$ss = gsub("\\.gz$", "", df$ss)
df$roi.nii = gsub("\\.gz$", "", df$roi.nii)

df$outfile = gsub("\\.nii", "_Masked.nii", df$raw)

save(df, file=file.path(rootdir, "Registration", 
  "All_Registration_Image_Names.Rda"))
all.df = df

ids.111 = read.csv(file.path(basedir, "111_patients.csv"), 
  stringsAsFactors= FALSE)
uid = ids.111$patientName
all.ids = ids.111$id

############################
# Take out 111 already ran for paper - this is for testing
df = df[ !df$id %in% all.ids, ]
rownames(df) =  NULL
############################

save(df, file=file.path(rootdir, "Registration", 
  "Registration_Image_Names.Rda"))

##########################################
# Copy files, unzip them, then mask the files
##########################################
melted = melt(xdf[, c("id", "raw", "copydir", "ss", "roi.nii")], 
  id.vars = c("id", "copydir"))
melted = melted[, c("value", "copydir")]

##########################################
## Copying files over
##########################################
overwrite = TRUE
m_ply(function(value, copydir) {
  del.file = file.path(copydir, 
    gsub(".nii.gz", ".nii", basename(value), fixed=TRUE)
  )
  if (file.exists(del.file)) file.remove(del.file)
  file.copy(value, copydir, overwrite=overwrite)
}, .data = melted, .progress = "text")

##########################################
## Unzipping files - for SPM
##########################################
melted$image = file.path(melted$copydir, basename(melted$value))
l_ply(.data=melted$image, .fun =function(x) {
  if (file.exists(x)){
    system(sprintf('gunzip --force "%s"', x))
  }
}, .progress = "text")
# }## end for loop

##########################################
## Making masked files
##########################################
Sys.setenv(FSLOUTPUTTYPE = "NIFTI")
m_ply(.data=all.df[, c("raw", "ss", "outfile")], 
  .fun = function(raw, ss, outfile) {
    if (!file.exists(outfile) | rerun){
      fslmask(file=raw, mask=ss, outfile=outfile)
    }
}, .progress = "text")



