#################################
# Regressions with % of ROI
# Author: John Muschelli
#################################
rm(list=ls())
library(cttools)
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


whichdir = "reoriented"
outcome = "NIHSS"

mses = NULL
outcomes = c("NIHSS", "GCS")
for (outcome in outcomes){
  load(file = file.path(basedir, 
    paste0("Cross_Validated_", outcome, "_Results.Rda"))
  )
  mses = cbind(mses, cm)
  rownames(mses) = names(cm)
}
colnames(mses) = outcomes

t(t(mses) < mses["Clot_Location_RC",])