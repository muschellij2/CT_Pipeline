###########################################
## This code creates a filename rda
##
## Author: John Muschelli
## Last updated: May 20, 2014
##############################################
##############################################
rm(list=ls())
library(plyr)
library(cttools)
library(devtools)
library(ROCR)
library(fslr)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

rerun = FALSE
outdir = file.path(basedir, "results")

iddirs = list.dirs(basedir, recursive=FALSE)
iddirs = iddirs[ grep("^\\d", basename(iddirs))]

files = unlist(sapply(iddirs, function(x){
  list.files(x, recursive = FALSE, 
    full.names=TRUE,
    pattern = ".*[.]nii")
  }))
names(files) = NULL
files = files[ !grepl("205-519_20110531_1140", 
  files)]
# files = list.files(basedir, recursive=TRUE, 
#   full.names=TRUE,
#            pattern=".*[.]nii")


fnames = basename(files)
fnames = sub("[.]nii", "ROI.nii", fnames)
# fnames = grep("^bws", fnames, value= TRUE)
# fnames = gsub("^bws", "", fnames)
# fnames = paste0(fnames, ".gz")
ids = gsub("(\\d\\d\\d-(\\d|)\\d\\d\\d)_.*", 
  "\\1", fnames)
fdf = data.frame(id = ids, 
  stringsAsFactors= FALSE)
fdf$iddir = file.path(basedir, fdf$id)
fdf$outdir = file.path(fdf$iddir, "Predictors")
makedir = sapply( fdf$outdir, dir.create, 
  showWarnings =FALSE)
fdf$roi = file.path(rootdir, "ROI_data", 
  fdf$id, fnames)
fdf$img = file.path(fdf$iddir, 
  gsub("ROI\\.nii", ".nii", 
  fnames))
fdf$mask = file.path(fdf$iddir, 
	"Skull_Stripped", 
	gsub("ROI\\.nii", "_SS_0.01_Mask.nii", 
    fnames))


fdf$ssimg = gsub("_Mask[.]nii[.]gz", 
  ".nii.gz", 
  fdf$mask)
fdf$n3dir = file.path(fdf$iddir, "N3Correct")

fdf$n3ssimg = gsub("[.]nii[.]gz", 
  "_N3Correct.nii.gz", 
    basename(fdf$ssimg))
fdf$n3ssimg = file.path(fdf$n3dir, fdf$n3ssimg)

fdf$n3img = gsub("[.]nii[.]gz", 
  "_N3Correct.nii.gz", 
    basename(fdf$img))
fdf$n3img = file.path(fdf$n3dir, fdf$n3img)

fdf$n3bias_field = gsub("_N3Correct", 
  "_BiasField", 
  fdf$n3img)
fdf$n3ssbias_field= gsub("_N3Correct", 
  "_BiasField", 
  fdf$n3ssimg)

fdf$n4img = gsub("_N3Correct[.]nii[.]gz", 
	"_N4Correct.nii.gz", 
    fdf$n3img)
fdf$n4ssimg = gsub("_N3Correct[.]nii[.]gz", 
	"_N4Correct.nii.gz", 
    fdf$n3ssimg)

fdf$n4bias_field = gsub("_N3Correct", 
  "_BiasField_N4", 
  fdf$n3img)
fdf$n4ssbias_field= gsub("_N3Correct", 
  "_BiasField_N4", 
  fdf$n3ssimg)


fdf$syndir = file.path(fdf$iddir, 
  "SyN_Registered")
fdf$synssimg = gsub("[.]nii[.]gz", 
  "_SyN.nii.gz", 
    basename(fdf$ssimg))
fdf$synssimg = file.path(fdf$syndir, 
  fdf$synssimg)
fdf$synssroi = gsub("_SyN[.]nii[.]gz", 
  "_SyN_ROI.nii.gz", 
  fdf$synssimg)
fdf$synssmask = gsub("_SyN[.]nii[.]gz", 
  "_SyN_Mask.nii.gz", 
  fdf$synssimg)

fdf$sinc_synssimg = gsub("_SyN[.]nii[.]gz", 
  "_SyN_sinc.nii.gz", 
  fdf$synssimg)
fdf$sinc_synssroi = gsub("_SyN[.]nii[.]gz", 
  "_SyN_sinc_ROI.nii.gz", 
  fdf$synssimg)
fdf$sinc_synssmask = gsub("_SyN[.]nii[.]gz", 
  "_SyN_sinc_Mask.nii.gz", 
  fdf$synssimg)


fdf$aff_dir = file.path(fdf$iddir, 
  "Affine_Registered")
fdf$aff_ssimg = gsub("_SyN[.]nii[.]gz", 
  "_Affine.nii.gz", 
  basename(fdf$synssimg))
fdf$aff_ssimg = file.path(fdf$aff_dir, 
  fdf$aff_ssimg)
fdf$aff_ssroi = gsub("_Affine[.]nii[.]gz", 
  "_Affine_ROI.nii.gz", 
  fdf$aff_ssimg)
fdf$aff_ssmask = gsub("_Affine[.]nii[.]gz", 
  "_Affine_Mask.nii.gz", 
  fdf$aff_ssimg)

fdf$sinc_aff_ssimg = gsub("_SyN[.]nii[.]gz", 
  "_Affine_sinc.nii.gz", 
  basename(fdf$synssimg))
fdf$sinc_aff_ssimg = file.path(fdf$aff_dir, 
  fdf$sinc_aff_ssimg)
fdf$sinc_aff_ssroi = gsub("_Affine[.]nii[.]gz", 
  "_Affine_sinc_ROI.nii.gz", 
  fdf$aff_ssimg)
fdf$sinc_aff_ssmask = gsub("_Affine[.]nii[.]gz", 
  "_Affine_sinc_Mask.nii.gz", 
  fdf$aff_ssimg)


fdf$rig_dir = file.path(fdf$iddir, 
  "Rigid_Registered")
fdf$rig_ssimg = gsub("_SyN[.]nii[.]gz", 
  "_Rigid.nii.gz", 
  basename(fdf$synssimg))
fdf$rig_ssimg = file.path(fdf$rig_dir, 
  fdf$rig_ssimg)

fdf$rig_ssroi = gsub("_Rigid[.]nii[.]gz", 
  "_Rigid_ROI.nii.gz", 
  fdf$rig_ssimg)
fdf$rig_ssmask = gsub("_Rigid[.]nii[.]gz", 
  "_Rigid_Mask.nii.gz", 
  fdf$rig_ssimg)

fdf$sinc_rig_ssimg = gsub("_SyN[.]nii[.]gz", 
  "_Rigid_sinc.nii.gz", 
  basename(fdf$synssimg))
fdf$sinc_rig_ssimg = file.path(fdf$rig_dir, 
  fdf$sinc_rig_ssimg)

fdf$sinc_rig_ssroi = gsub(
  "_Rigid_sinc[.]nii[.]gz", 
  "_Rigid_sinc_ROI.nii.gz", 
  fdf$sinc_rig_ssimg)
fdf$sinc_rig_ssmask = gsub(
  "_Rigid_sinc[.]nii[.]gz", 
  "_Rigid_sinc_Mask.nii.gz", 
  fdf$sinc_rig_ssimg)


############################
# Only 1 scan right now
############################
fdf = ddply(fdf, .(id), function(x){
  x[1,]
})

ids = c("100-318", "100-362", "100-365", 
  "101-306", 
  "101-307", "101-308", 
"102-317", "102-322", "102-323", "102-324", 
"102-331", 
"102-360", 
"102-367", "102-374", "102-391", "102-393", 
"102-403", 
"102-406", 
"120-376", "131-310", "131-316", "131-334", 
"131-354", 
"133-409", 
"133-417", "134-304", "134-305", "134-320", 
"134-327", 
"134-345", 
"134-380", "134-381", "134-382", "134-392", 
"134-408", 
"134-412", 
"134-416", "152-302", "152-303", "152-353", 
"157-328", 
"157-329", 
"157-332", "157-335", "157-336", "157-370", 
"157-372", 
"157-399", 
"157-410", "161-413", "173-312", "173-313", 
"173-325", 
"173-341", 
"173-361", "173-364", "173-368", "173-384",
 "173-396", 
"173-404", 
"175-387", "175-397", "175-405", "179-373",
 "179-383", 
"179-386", 
"179-394", "179-395", "179-402", "184-388", 
"191-301", 
"191-311", 
"191-314", "191-315", "191-319", "191-321", 
"191-333", 
"191-375", 
"191-400", "205-509", "205-517", "205-519", 
"216-390", 
"219-350", 
"222-337", "222-357", "222-358", "223-355", 
"223-369", 
"223-407", 
"225-502", "225-503", "225-504", "225-505", 
"225-506", 
"225-507", 
"225-510", "225-511", "225-515", "225-524", 
"230-356", 
"230-363", 
"230-366", "230-371", "230-377", "232-514", 
"232-516", 
"234-385", 
"265-389", "265-398", "289-518", "289-525")

fdf = fdf[ fdf$id %in% ids, ]




###########################################
# Create the groupings
###########################################
set.seed(20150504)
non.aggmods = 10
vec = seq(length(ids))
ind = sample(1:length(vec), size=non.aggmods)
train.ind = sort(vec[ind])
vec = vec[-ind]
nvalid = ceiling(length(vec)/2)
ind = sample(1:length(vec), size=nvalid)
valid.ind = sort(vec[ind])
test.ind = sort(vec[-ind])

df = data.frame(id = ids)
df$group = NA
df$group[train.ind] = "Train"
df$group[valid.ind] = "Validation"
df$group[test.ind] = "Test"

fdf = merge(fdf, df, by="id", sort=FALSE)
if (any(is.na(fdf$group))){
  stop(paste0("Something went wrong with ", 
    "group designation"))
}

stopifnot(all(file.exists(
  unlist(fdf[, c("roi", "img", "ssimg")])))
)

outfile = file.path(outdir, 
  "111_Filenames.Rda")
save(fdf, fnames, file = outfile)


cn = c("id", "iddir", "outdir", "roi", 
  "img", "mask", "ssimg", "rig_dir", 
  "rig_ssimg", 
"rig_ssroi", "rig_ssmask",
"group")
fdf = fdf[, cn]
fdf$preddir = fdf$outdir
outfile = file.path(outdir, 
  "Reseg_111_Filenames.Rda")
save(fdf, fnames, file = outfile)
