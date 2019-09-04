####################################################################
## This code creates a filename rda for all images
##
## Author: John Muschelli
## Last updated: Dec 7, 2014
####################################################################
####################################################################
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

files = list.files(basedir, recursive=TRUE, full.names=TRUE,
           pattern=".*[.]nii")
fnames = basename(files)
fnames = grep("^bws", fnames, value= TRUE)
fnames = gsub("^bws", "", fnames)
fnames = paste0(fnames, ".gz")
ids = gsub("(\\d\\d\\d-(\\d|)\\d\\d\\d)_.*", "\\1", fnames)
fdf = data.frame(id = ids, stringsAsFactors= FALSE)
fdf$iddir = file.path(basedir, fdf$id)
fdf$outdir = file.path(fdf$iddir, "Predictors")
makedir = sapply( fdf$outdir, dir.create, showWarnings =FALSE)
fdf$roi = file.path(rootdir, "ROI_data", fdf$id, fnames)
fdf$img = file.path(fdf$iddir, gsub("ROI\\.nii", ".nii", fnames))
fdf$mask = file.path(fdf$iddir, 
	"Skull_Stripped", 
	gsub("ROI\\.nii", "_SS_0.01_Mask.nii", fnames))


fdf$ssimg = gsub("_Mask[.]nii[.]gz", ".nii.gz", fdf$mask)
fdf$n3dir = file.path(fdf$iddir, "N3Correct")

fdf$n3ssimg = gsub("[.]nii[.]gz", "_N3Correct.nii.gz", 
    basename(fdf$ssimg))
fdf$n3ssimg = file.path(fdf$n3dir, fdf$n3ssimg)

fdf$n3img = gsub("[.]nii[.]gz", "_N3Correct.nii.gz", 
    basename(fdf$img))
fdf$n3img = file.path(fdf$n3dir, fdf$n3img)

fdf$n3bias_field = gsub("_N3Correct", "_BiasField", fdf$n3img)
fdf$n3ssbias_field= gsub("_N3Correct", "_BiasField", fdf$n3ssimg)

fdf$n4img = gsub("_N3Correct[.]nii[.]gz", 
	"_N4Correct.nii.gz", 
    fdf$n3img)
fdf$n4ssimg = gsub("_N3Correct[.]nii[.]gz", 
	"_N4Correct.nii.gz", 
    fdf$n3ssimg)

fdf$n4bias_field = gsub("_N3Correct", "_BiasField_N4", 
  fdf$n3img)
fdf$n4ssbias_field= gsub("_N3Correct", "_BiasField_N4", 
  fdf$n3ssimg)


fdf$syndir = file.path(fdf$iddir, "SyN_Registered")
fdf$synssimg = gsub("[.]nii[.]gz", "_SyN.nii.gz", 
    basename(fdf$ssimg))
fdf$synssimg = file.path(fdf$syndir, fdf$synssimg)
fdf$synssroi = gsub("_SyN[.]nii[.]gz", "_SyN_ROI.nii.gz", 
  fdf$synssimg)
fdf$synssmask = gsub("_SyN[.]nii[.]gz", "_SyN_Mask.nii.gz", 
  fdf$synssimg)

fdf$sinc_synssimg = gsub("_SyN[.]nii[.]gz", "_SyN_sinc.nii.gz", 
  fdf$synssimg)
fdf$sinc_synssroi = gsub("_SyN[.]nii[.]gz", 
  "_SyN_sinc_ROI.nii.gz", 
  fdf$synssimg)
fdf$sinc_synssmask = gsub("_SyN[.]nii[.]gz", 
  "_SyN_sinc_Mask.nii.gz", 
  fdf$synssimg)


fdf$aff_dir = file.path(fdf$iddir, "Affine_Registered")
fdf$aff_ssimg = gsub("_SyN[.]nii[.]gz", "_Affine.nii.gz", 
  basename(fdf$synssimg))
fdf$aff_ssimg = file.path(fdf$aff_dir, fdf$aff_ssimg)
fdf$aff_ssroi = gsub("_Affine[.]nii[.]gz", "_Affine_ROI.nii.gz", 
  fdf$aff_ssimg)
fdf$aff_ssmask = gsub("_Affine[.]nii[.]gz", "_Affine_Mask.nii.gz", 
  fdf$aff_ssimg)

fdf$sinc_aff_ssimg = gsub("_SyN[.]nii[.]gz", 
  "_Affine_sinc.nii.gz", 
  basename(fdf$synssimg))
fdf$sinc_aff_ssimg = file.path(fdf$aff_dir, fdf$sinc_aff_ssimg)
fdf$sinc_aff_ssroi = gsub("_Affine[.]nii[.]gz", 
  "_Affine_sinc_ROI.nii.gz", 
  fdf$aff_ssimg)
fdf$sinc_aff_ssmask = gsub("_Affine[.]nii[.]gz", 
  "_Affine_sinc_Mask.nii.gz", 
  fdf$aff_ssimg)


fdf$rig_dir = file.path(fdf$iddir, "Rigid_Registered")
fdf$rig_ssimg = gsub("_SyN[.]nii[.]gz", "_Rigid.nii.gz", 
  basename(fdf$synssimg))
fdf$rig_ssimg = file.path(fdf$rig_dir, fdf$rig_ssimg)

fdf$rig_ssroi = gsub("_Rigid[.]nii[.]gz", "_Rigid_ROI.nii.gz", 
  fdf$rig_ssimg)
fdf$rig_ssmask = gsub("_Rigid[.]nii[.]gz", "_Rigid_Mask.nii.gz", 
  fdf$rig_ssimg)

fdf$sinc_rig_ssimg = gsub("_SyN[.]nii[.]gz", 
  "_Rigid_sinc.nii.gz", 
  basename(fdf$synssimg))
fdf$sinc_rig_ssimg = file.path(fdf$rig_dir, fdf$sinc_rig_ssimg)

fdf$sinc_rig_ssroi = gsub("_Rigid_sinc[.]nii[.]gz", 
  "_Rigid_sinc_ROI.nii.gz", 
  fdf$sinc_rig_ssimg)
fdf$sinc_rig_ssmask = gsub("_Rigid_sinc[.]nii[.]gz", 
  "_Rigid_sinc_Mask.nii.gz", 
  fdf$sinc_rig_ssimg)


############################
# Only 1 scan right now
############################
fdf = ddply(fdf, .(id), function(x){
  x[1,]
})

outfile = file.path(outdir, "111_Filenames.Rda")
save(fdf, fnames, file = outfile)
