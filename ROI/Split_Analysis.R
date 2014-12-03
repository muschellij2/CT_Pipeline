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
                                levels= 
                                  c("Lobar", "Globus Pallidus", 
                                    "Putamen", "Thalamus"))
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
# Over 10 needs to be in sample - not outcome dep 
###############################################
ncut = 10
dim(mat)
sum(rs > ncut)
mat = mat[rs > ncut, ]

###############################################
# Cross validation folds
###############################################
# set.seed(20141106)
N <- nrow(demog)
nfolds <- 2
folds <-c(sapply(1:nfolds, rep, length=ceiling(N/nfolds)))
folds <- folds[1:N]
B = 10
l.pval = l.adjr2 = vector(mode="list", length=B)
for (iB in 1:B){
  folds = sample(folds)
  
  ###############################################
  # Cross validation Over 10 needs to be in sample
  ###############################################
  top.ords = c(1000, 2000, 3000)
  pvals = c(0.05, .01, .001)
  
  runnames = c(paste0("pval.", pvals), 
               paste0("ord.", top.ords))
  
  ifold = 1
  res.adjs = res.mats = NULL
  for (ifold in sort(unique(folds))){
    fold.ind = which(folds == ifold)
    outfold.ind = which(folds != ifold)
    
    ###############################################
    # Get corresponding binary image values
    ###############################################
    runmat = t(mat[,outfold.ind ])
    class(runmat) = "numeric"
    testmat = t(mat[,fold.ind ])
    class(testmat) = "numeric"  
    
    ###############################################
    # Get Y values
    ###############################################  
    test.demog = demog[fold.ind,]
    
    run.demog = demog[outfold.ind,]
    runY = run.demog$Y
    mytime = system.time({
      mods = fast_lm(runY, X=runmat, Z = NULL, 
                     spot.check = TRUE, ncheck = 10, verbose=FALSE)
    })  
    mods$AIC = extractAIC(mods)
    
    ###############################################
    # P-value Cutoff
    ###############################################  
    for (pval in pvals){
      ind = which(mods$p.val <= pval)
      rr = rowMeans(runmat[, ind, drop=FALSE])
      if (all(is.nan(rr))){
        rr = 0
      }
      run.demog[, paste0("pval.", pval)] = rr
      
      rr = rowMeans(testmat[, ind, drop=FALSE])
      if (all(is.nan(rr))){
        rr = 0
      }
      test.demog[, paste0("pval.", pval)] = rr
    }
    
    ###############################################
    # Top P-value Cutoff
    ###############################################    
    for (top.ord in top.ords){
      ord = order(mods$p.val)
      ind = which(ord <= top.ord)
      rr = rowMeans(runmat[, ind, drop=FALSE])
      if (all(is.nan(rr))){
        rr = 0
      }      
      run.demog[, paste0("ord.", top.ord)] = rr
      
      rr = rowMeans(testmat[, ind, drop=FALSE])
      if (all(is.nan(rr))){
        rr = 0
      }      
      test.demog[, paste0("ord.", top.ord)] = rr      
    }
    
    
    ###############################################
    # Likelihood ratio tests for models with the location in there and HPR
    ###############################################     
    loc.mod = lm(Y ~ ., 
                 data= test.demog[, 
                                  c("Y", "Age", "Sex", "Diagnostic_ICH", 
                                    "Clot_Location_RC")])      
    
    null.mod = lm(Y ~ ., 
                 data= test.demog[, 
                                  c("Y", "Age", "Sex", "Diagnostic_ICH")])  
    iname = runnames[1]
    res.adj2 = rep(NA, length=length(runnames)+1)
    names(res.adj2) = c(runnames, "Clot_Location_RC")
    res.mat = matrix(NA, nrow=3, ncol=length(runnames))
    colnames(res.mat) = runnames
    rownames(res.mat) = c("Location LR Test P-value", "HPR LR Test P-value", 
                          "HPR LR Test Null P-value")
    res.adj2['Clot_Location_RC'] = summary(loc.mod)$adj.r.squared
    for (iname in runnames){
      full.mod = lm(Y ~ ., 
                    data= test.demog[, 
                                     c("Y", "Age", "Sex", "Diagnostic_ICH", 
                                       "Clot_Location_RC", iname)])
      pval.mod = lm(Y ~ ., 
                    data= test.demog[, 
                                     c("Y", "Age", "Sex", "Diagnostic_ICH", 
                                       iname)])    
      loc.fpval = anova(full.mod, pval.mod)$"Pr(>F)"[2]
      pval.fpval = anova(full.mod, loc.mod)$"Pr(>F)"[2]
      pval.null = anova(null.mod, pval.mod)$"Pr(>F)"[2]
      res.adj2[iname] = summary(pval.mod)$adj.r.squared
      
      res.mat["Location LR Test P-value", iname] = loc.fpval
      res.mat["HPR LR Test P-value", iname] = pval.fpval
      res.mat["HPR LR Test Null P-value", iname] = pval.null
      
    }
    rownames(res.mat) = paste0(rownames(res.mat), ".", ifold)
    res.mats = rbind(res.mats, res.mat)
    res.adjs = rbind(res.adjs, res.adj2)
#     print(ifold)
  }
  
  rownames(res.adjs) = sort(unique(folds))
  
#   print(res.mats)
#   print(res.adjs)
  l.pval[[iB]] = res.mats
  l.adjr2[[iB]] = res.adjs


}



print(l.pval)