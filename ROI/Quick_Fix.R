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
demog$Clot_Location_RC = factor(demog$Clot_Location_RC, 
  levels= c("Lobar", "Globus Pallidus", "Putamen", "Thalamus"))
demog$LOC = demog$Clot_Location_RC

zform = ~ Age + Sex + Diagnostic_ICH
Z = model.matrix(object = zform, data = demog)
cc = complete.cases(Z) & complete.cases(demog$Y)
demog = demog[cc,]
Z = model.matrix(object = zform, data = demog)

outfile = file.path(outdir, "Voxel_Matrix.Rda")
load(file=outfile )


B = 1000
Y = sapply(seq(B), function(x){
  sample(demog$Y)
})


#### keeping if over 10 people have ICH in that locaiton
ncut = 10
all.nvox = sum(rs > 0)

dim(mat)
sum(rs > ncut)
runmat = mat[rs > ncut, ]

orig.X = t(runmat)
class(orig.X) = "numeric"

X = t(runmat[, cc, drop=FALSE])
class(X) = "numeric"

cn = colnames(Z)
ZZ = Z[, !cn %in% c("(Intercept)")]
z = lm(Y ~ X[,5] + ZZ)
z.noint = lm(Y ~ X[,5] + Z - 1)


verbose= TRUE
tol = 1e-12
ncheck = 10

mytime = system.time({
  mods = fast_lm(Y, X, Z, spot.check = TRUE, ncheck = 10)
})


pval = .01
top.ord = 3000

under.pval = mods$p.val <= pval
l.pval = apply(under.pval, 2, which)

l.ord = alply(mods$p.val, .margin=2, function(x){
  ord = order(x)
  which(ord <= top.ord)
}, .progress= "text")


# cov.pval = sapply(l.pval, function(x){
#   colMeans(mat[x, ])
# })

cov.ord = sapply(l.ord, function(x){
  colMeans(runmat[x, ])
})
# pvals 

new.cov = cov.ord[cc, ]
orig.Y = demog$Y
permute.mod =  fast_lm(orig.Y, X=new.cov, Z, spot.check = TRUE, 
  ncheck = 10)

################
# Check against truth

x = load(file.path(rootdir, "CT_Pipeline/Top_3000_Pvalues_df.Rda"))

##### subsetting the correct patients and then getting means
submat = submat[,cc] 
perc = colMeans(submat)

#### should double check against IDs
demog$perc_ROI = perc * 100

#############
# Making the formula for the adjusted analysis
#############
form = as.character(zform)
form = c(paste0("Y", form[1]), paste0(form[2], "+ perc_ROI"))
form = as.formula(form)
reg.mod = lm(formula=form, data=demog)

#################
# Checking the p-values here - univarite p-values
#################
check.mod = fast_lm(orig.Y, X= orig.X, Z = NULL, spot.check=TRUE)

##### getting top values
ord = order(check.mod$p.val)
rrn = which(rs > ncut)
rrn = rrn[ord]
myrn = rrn[seq(top.ord)]

ssmat = mat[ rn, cc]
pp = colMeans(ssmat)




real.mod = fast_lm(orig.Y, X= as.matrix(perc), Z, spot.check=TRUE)
real.r2 = real.mod$r.squared
d.check = lm(orig.Y ~ perc + Z)
check.r2 = summary(d.check)$r.squared
hist(permute.mod$r.squared)
