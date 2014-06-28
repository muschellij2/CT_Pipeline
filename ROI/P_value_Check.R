#################################
# Regressions with % of ROI
# Author: John Muschelli
#################################
rm(list=ls())
library(oro.nifti)
library(plyr)
library(scales)
library(RColorBrewer)
library(data.table)
library(cttools)
library(fslr)
library(ggplot2)
library(grid)
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


whichdir = "reoriented"
outcome = "NIHSS"
adder = paste0(outcome, "_")
if (outcome == "NIHSS"){
  adder = ""
}

get.id = function(x){
  ss = strsplit(x, "_")
  ss = sapply(ss, head, 1)
  ss = gsub(".*(\\d\\d\\d-.*)", "\\1", ss)
  ss
}

id_to_pname = function(x){
  as.numeric(gsub("-", "", x))
}

nkeeps = c(1000, 2000, 3000, .001, 0.01, 0.05)
cn = c("With_Perc", "With_Clot", "Null", 
       "Nosex_With_Perc", "Nosex_With_Clot", 
       "Nosex_Null", "nkeep", "pval", "N_Gr0")


demog = read.csv(file=file.path(basedir, "Demog_NIHSS_Mask.csv"), 
                 stringsAsFactors=FALSE)
demog$Base_ICH_10 = demog$Diagnostic_ICH /10

measures  = c("adj.r.squared", "r.squared", "sigma")
reses = vector("list", length=length(measures))
names(reses) = measures
epics = aics = reses 

if (outcome == "GCS") {
  demog$Y = demog$Enrollment_GCS_Add
} else if (outcome == "NIHSS"){
  demog$Y = demog$Enrollment_NIHSS_Total
} else {
  stop(paste0("Outcome ", outcome, " not implemented"))
}


demog$Clot_Location_RC = gsub("Palidus", "Pallidus", 
  demog$Clot_Location_RC )
demog$Clot_Location_RC = factor(demog$Clot_Location_RC, levels= c("Lobar", "Globus Pallidus", "Putamen", "Thalamus"))
demog$LOC = demog$Clot_Location_RC

outfile = file.path(outdir, "Voxel_Matrix.Rda")
load(file=outfile )

#### keeping if over 10 people have ICH in that locaiton
ncut = 10
all.nvox = sum(rs > 0)

dim(mat)
sum(rs > ncut)
mat = mat[rs > ncut, ]


X = t(mat[, drop=FALSE])
class(X) = "numeric"


# test = X[, 1:100]
test = X
imod = 1
for (imod in seq(ncol(test))){
  x = test[, imod]
  m1 = lm(Y ~ Age + Sex + Diagnostic_ICH + x, data=demog) 
  pv1 =  coef(summary(m1))["x",  "Pr(>|t|)"]
  m2 = lm(x ~ Age + Sex + Diagnostic_ICH + Y, data=demog) 
  pv2 =  coef(summary(m2))["Y",  "Pr(>|t|)"]  
  stopifnot(abs(pv1 - pv2) < 1e-12)
  print(imod)
}
