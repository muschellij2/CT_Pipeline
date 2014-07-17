#################################
# Regressions with % of ROI
# Author: John Muschelli
#################################
rm(list=ls())
library(cttools)
library(oro.nifti)
library(scales)
library(RColorBrewer)
library(data.table)
library(ggplot2)
library(grid)
library(plyr)
homedir = "/Applications"
rootdir = "~/CT_Registration"
basedir = file.path(rootdir, "data")
outdir = basedir
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
  basedir = file.path(rootdir, "Registration")
  outdir = file.path(basedir, "results")
}
progdir = file.path(rootdir, "programs")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")


get.id = function(x){
  ss = strsplit(x, "_")
  ss = sapply(ss, head, 1)
  ss = gsub(".*(\\d\\d\\d-.*)", "\\1", ss)
  ss
}

id_to_pname = function(x){
  as.numeric(gsub("-", "", x))
}

demog = read.csv(file=file.path(basedir, "Demog_NIHSS_Mask.csv"), 
 stringsAsFactors=FALSE)
demog$Base_ICH_10 = demog$Diagnostic_ICH /10


template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)
dtemp = dim(temp)

#### load voxel data
outfile = file.path(outdir, "Voxel_Matrix.Rda")
vmat = load(file=outfile )

stopifnot(nrow(mat) == prod(dtemp))

nim = temp
nim[is.na(nim)] = 0
nim[!is.na(nim)] = 1

max.slice = dtemp[1] 
mid.slice = (max.slice+1)/2

w = which(nim > 0, arr.ind=TRUE)
w2 = which(nim > 0)
w = cbind(w, ind=w2)
w = data.frame(w)
w = w[ w$dim1 != mid.slice, ]
### remember left / right for oro.nifti
w$left = w$dim1 > mid.slice

i.left = w$ind[w$left]
i.right = w$ind[!w$left]

s.left = colSums(mat[i.left, ])
s.right = colSums(mat[i.right, ])

stopifnot(!any(s.left == s.right))

df = data.frame(cbind(left=s.left, right=s.right))
df$patientName = id_to_pname(get.id(names(s.left)))
df$side = ifelse(df$left > df$right, "Left", "Right")

mm = merge(demog[, c("patientName", "Hemisphere")], df, all=TRUE)

table(mm$side, mm$Hemisphere)
