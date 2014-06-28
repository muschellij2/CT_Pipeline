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


demog$Clot_Location_RC = gsub("Palidus", "Pallidus", demog$Clot_Location_RC )
demog$Clot_Location_RC = factor(demog$Clot_Location_RC, levels= c("Lobar", "Globus Pallidus", "Putamen", "Thalamus"))
demog$LOC = demog$Clot_Location_RC

Z = model.matrix(~ Age + Sex + Diagnostic_ICH, data= demog)
cc = complete.cases(Z) & complete.cases(demog$Y)
demog = demog[cc,]
Z = model.matrix(~ Age + Sex + Diagnostic_ICH, data= demog)

B = 10
Y = sapply(seq(B), function(x){
  sample(demog$Y)
})

outfile = file.path(outdir, "Voxel_Matrix.Rda")
load(file=outfile )

#### keeping if over 10 people have ICH in that locaiton
ncut = 10
all.nvox = sum(rs > 0)

dim(mat)
sum(rs > ncut)
mat = mat[rs > ncut, ]


X = t(mat[, cc, drop=FALSE])
class(X) = "numeric"

# Y = Y[cc,, drop=FALSE]
# Z = Z[cc,, drop=FALSE]
n = nrow(Y)

#################
# Get residualls on Z
#################
ZtZ = crossprod(Z)
ZtZ.inv = solve(ZtZ)
Z.H = diag(n) - Z %*% tcrossprod(ZtZ.inv, Z)

# Q = aaply(X, 2, function(x){
#   r = crossprod(x, Z.H)
#   xx = 1/(r %*% x)
#   xx = xx %*% r
#   as.vector(xx)
# }, .progress="text")

# betas = Q %*% Y

X = Z.H %*% X
Y = Z.H %*% Y

rxvars = colVars(X)
rxss = colSums(X^2)

cov_XY = crossprod(X, Y)
beta_x = cov_XY / rxss

ses = matrix(NA, 
    nrow=nrow(beta_x), 
    ncol=B)
stopifnot(nrow(X) == n)
stopifnot(nrow(beta_x)==  ncol(X))
iB = 2
for (iB in seq(B)){
  b = beta_x[, iB]
  system.time({
    m = matrix(b, ncol=ncol(X), nrow=n, byrow=TRUE)
    yhat = m * X
  })
  # ###### multiplication goes column-wise 
  # system.time({
  #   xx = X * b
  # })


  check = X[,1] * b[1]
  stopifnot(all(abs(yhat[,1] - check) == 0))
  # check[,1] - xx[,1]

  r = yhat - Y[,iB]
  ### double check worked correctly
  check = (yhat[,1] - Y[,iB])
  stopifnot(all(abs(r[,1] - check) == 0))


  rss= colSums( ( r ) ^2)
  p = 5
  bse = sqrt(rxss * rss/ (n-p))
  ses[, iB] = bse
  print(iB)
}

xx = X[,1]
m1 = lm(Y[, 1] ~ Z + xx - 1)
bcoef = coef(m1)
abs(bcoef["xx"] - beta_x[1, 1]) < 1e-12
smod = summary(m1)
se = coef(smod)["xx", "Std. Error"]
ses[1,1]
# iB = 1
# ses = matrix(NA, 
#     nrow=nrow(beta_x), 
#     ncol=ncol(beta_x))
# for (iB in seq(B)){
#   b = beta_x[,iB]
#   yhat = t(b * X)
#   ses[,iB] = sqrt(rowSums( (Y[, iB] - yhat)^2 / (n- p)))
#   print(iB)
# }



iX = 1
iY = 1
xx = X[,iX]
ZZ = Z[, !(colnames(Z) %in% "(Intercept)")]
check = lm(Y[,iY] ~ ZZ + xx)
d = abs(coef(check)[["xx"]] - beta_x[iX, iY])
stopifnot(d < 1e-12)
checkse = coef(summary(check))["xx", "Std. Error"]
bse[,1]

# cov_RXY = crossprod(R_X, R_Y)
# beta_x = cov_RXY / rxss
# y_hats = R_X %*% beta_x
# cs = aaply(beta_x, .margins=2, function(b){
#   xx = (b * R_X)^2
#   colSums(xx)
# }, .progress = "text")
# resid.var = colSums(y_hats^2)



# for (i in 1:1000){
#   iX = sample(seq(ncol(X)), size=1)
#   iY = sample(seq(ncol(Y)), size=1)
#   xx = X[,iX]
#   check = lm(Y[,iY] ~ ZZ + xx)
#   d = abs(coef(check)[["xx"]] - beta_x[iX, iY])
#   stopifnot(d < 1e-12)
# }
#     nkeep = nkeeps[ikeep]
