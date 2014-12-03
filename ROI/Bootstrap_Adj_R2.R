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
# outcome = "GCS"
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

demog = read.csv(file=file.path(basedir, "Demog_NIHSS_Mask.csv"), 
                 stringsAsFactors=FALSE)
demog$Base_ICH_10 = demog$Diagnostic_ICH /10

if (outcome == "GCS") {
  demog$Y = demog$Enrollment_GCS_Add
} else if (outcome == "NIHSS"){
  demog$Y = demog$Enrollment_NIHSS_Total
} else {
  stop(paste0("Outcome ", outcome, " not implemented"))
}


demog$Clot_Location_RC = gsub("Palidus", "Pallidus", 
                              demog$Clot_Location_RC )
demog$Clot_Location_RC = factor(demog$Clot_Location_RC, 
                                levels= c("Lobar", "Globus Pallidus", "Putamen", "Thalamus"))
demog$LOC = demog$Clot_Location_RC

cc = complete.cases(demog$Y)
demog = demog[cc,]


zform = ~ Age + Sex + Diagnostic_ICH
# Z = model.matrix(object = zform, data = demog)
# Z = model.matrix(object = zform, data = demog)

###############################################
# Load and subset matrix
###############################################
outfile = file.path(outdir, "Voxel_Matrix.Rda")
load( file=outfile )
mat = mat[, which(cc), drop=FALSE]

###############################################
# Cross validation folds
###############################################
set.seed(20141106)
N <- nrow(demog)
B = 1000
folds <-sapply(1:B, function(x) sample(N))

###############################################
# Over 10 needs to be in sample - not outcome dep 
###############################################
ncut = 10
dim(mat)
sum(rs > ncut)
mat = mat[rs > ncut, ]

###############################################
# Cross validation Over 10 needs to be in sample
###############################################
top.ords = c(1000, 2000, 3000)
pvals = c(0.05, .01, .001)

runnames = c(paste0("pval.", pvals), 
             paste0("ord.", top.ords))
r2.df = matrix(NA, ncol= length(runnames)+1, nrow=B)
colnames(r2.df)= c(runnames, "Clot_Location_RC")
r2.df = as.data.frame(r2.df)
iB = 1
for (iB in seq(B)){
  outfold.ind = folds[, iB]
  
  ###############################################
  # Get corresponding binary image values
  ###############################################
  runmat = t(mat[,outfold.ind ])
  class(runmat) = "numeric"

  
  ###############################################
  # Get Y values
  ###############################################    
  run.demog = demog[outfold.ind,]
  runY = run.demog$Y
  mytime = system.time({
    mods = fast_lm(runY, X=runmat, Z = NULL, 
                   spot.check = TRUE, ncheck = 10)
  })  
  
  ###############################################
  # P-value Cutoff
  ###############################################  
  for (pval in pvals){
    ind = which(mods$p.val <= pval)
    rr = rowMeans(runmat[, ind, drop=FALSE])
    run.demog[, paste0("pval.", pval)] = rr 
  }
  
  ###############################################
  # Top P-value Cutoff
  ###############################################    
  for (top.ord in top.ords){
    ord = order(mods$p.val)
    ind = which(ord <= top.ord)
    rr = rowMeans(runmat[, ind, drop=FALSE])
    run.demog[, paste0("ord.", top.ord)] = rr  
  }
  
  locform =  as.formula(c("Y ~",
                          paste0(
                            as.character(zform)[2], "+ Clot_Location_RC")
  )
  )
  
  runforms =  sapply(runnames, function(x) {
    as.formula(c("Y ~",
                 paste0(
                   as.character(zform)[2], "+ ", x)))
  })
  
  forms = c(locform, runforms)
  names(forms)[[1]] = "Clot_Location_RC"
  ########################
  # Run models
  ########################  
  preds = sapply(forms, function(formula){
    mod = lm(formula=formula, data=run.demog)
    smod = summary(mod)
    smod$adj.r.squared
  })
  
  cn = names(preds)
  r2.df[iB, cn] = preds
  print(iB)
}
