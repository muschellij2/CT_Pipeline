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
outcome = "GCS"
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

cc = complete.cases(demog$Y)
demog = demog[cc,]


zform = ~ Age + Sex + Diagnostic_ICH
# Z = model.matrix(object = zform, data = demog)
# Z = model.matrix(object = zform, data = demog)

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

X = t(runmat[, cc, drop=FALSE])
class(X) = "numeric"

verbose= TRUE
tol = 1e-12
ncheck = 10
Z = NULL

mytime = system.time({
  # mods = fast_lm(Y, X, Z, spot.check = TRUE, ncheck = 10)
  mods = fast_lm(Y, X, Z = NULL, spot.check = TRUE, ncheck = 10)
})


pval = .01
top.ords = c(1000, 2000, 3000)
pvals = c(0.05, .01, .001)

under.pval = lapply(pvals, function(pval) {
  x = mods$p.val <= pval
  l.x = apply(x, 2, which)
  l.len  = sapply(l.x, length)
  print(sum(l.len == 0))
  l.x
})

l.ord = llply(top.ords, function(top.ord) {
  alply(mods$p.val, .margin=2, function(x){
    ord = order(x)
    l = which(ord <= top.ord)
    # print(length(l))
    l
}, .progress= "text")
})

l.cov = laply(l.ord, function(ord) {
  laply(ord, function(x){
    rr = runmat[x, cc]
    if ( nrow(rr) == 0 ) {
      mns = rep(1, length=sum(cc))
      names(mns) = colnames(runmat)[cc]
      return(mns)
    }
    colMeans(rr)
}, .progress= "text")
})


cov.ord = sapply(l.ord, function(x){
  colMeans(runmat[x, ])
})

new.cov = cov.ord[cc, ]
orig.Y = demog$Y
permute.mod =  fast_lm(orig.Y, X=new.cov, Z=zform, data=demog,
  spot.check = TRUE, 
  ncheck = 10)

################
# Check against truth
###############


nkeeps = c(1000, 2000, 3000)

for (pval in c(0.05, .01, .001)){

  nkeeps = c(nkeeps, sum(yord <= pval))
}


for (nkeep in nkeeps){
  rn = rrn[seq(nkeep)]

  outstub = file.path(outdir, 
    paste0(adder, "Top_", nkeep, "_pvalues"))

  }
}



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

stopifnot(abs(pp-perc) < 1e-12)
stopifnot(all(submat == ssmat))


real.mod = fast_lm(orig.Y, X= as.matrix(perc), Z = zform, 
  data=demog, spot.check=TRUE)
real.r2 = real.mod$r.squared
check.r2 = summary(reg.mod)$r.squared
r.r2 = range(permute.mod$r.squared)
xlim = c(min(check.r2, r.r2), max(check.r2, r.r2))
hist(permute.mod$r.squared, xlim=xlim)
abline(v= check.r2)
