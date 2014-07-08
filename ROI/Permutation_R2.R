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
load( file=outfile )


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


########################
# Run models
########################

verbose= TRUE
tol = 1e-12
ncheck = 10
Z = NULL

mytime = system.time({
  # mods = fast_lm(Y, X, Z, spot.check = TRUE, ncheck = 10)
  mods = fast_lm(Y, X, Z = NULL, spot.check = TRUE, ncheck = 10)
})


top.ords = c(1000, 2000, 3000)
pvals = c(0.05, .01, .001)


under.pval = llply(pvals, function(pval) {
    x = mods$p.val <= pval
    l.x = apply(x, 2, which)
    l.len  = sapply(l.x, length)
    print(sum(l.len == 0))
    # print(length(l))
    l.x
}, .progress= "text")

names(under.pval) = pvals

l.pcov = llply(under.pval, function(ord) {
  laply(ord, function(x){
    rr = runmat[x, cc, drop=FALSE]
    if ( nrow(rr) == 0 ) {
      mns = rep(1, length=sum(cc))
      names(mns) = colnames(runmat)[cc]
      return(mns)
    }
    colMeans(rr)
  }, .progress= "text")
})

names(l.pcov) = pvals



l.ord = llply(top.ords, function(top.ord) {
  alply(mods$p.val, .margin=2, function(x){
    ord = order(x)
    l = which(ord <= top.ord)
    # print(length(l))
    l
}, .progress= "text")
})
names(l.ord) = top.ords

l.cov = llply(l.ord, function(ord) {
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

names(l.cov) = top.ords


all.l = c(l.pcov, l.cov)
# cov.ord = sapply(l.ord, function(x){
#   colMeans(runmat[x, ])
# })

all.run = c(pvals, top.ords)


orig.Y = demog$Y


#################
# Checking the p-values here - univarite p-values
#################

pv = load(file.path(outdir, paste0(adder, "Pvalue_Matrix.Rda")))
true.beta = res[,"Xbeta","mod.1"]
true.pval = res[,"X","mod.1"]

check.mod = fast_lm(orig.Y, X= X, Z = NULL, 
  spot.check=TRUE, 
  ncheck=100)

def.check = abs(check.mod$p.val - true.pval)
stopifnot(all(def.check < 1e-12))
def.check = abs(check.mod$coef - true.beta)
stopifnot(all(def.check < 1e-12))


irun = 1

pdfname = file.path(outdir, 
  paste0(adder, "Permutation_Distribution.pdf"))
pdf(pdfname)
for (irun in seq_along(all.run)){
  runmod = as.character(all.run[irun])

  # new.cov = cov.ord[cc, ]
  new.cov = t(all.l[[runmod]])
  permute.mod =  fast_lm(orig.Y, X=new.cov, Z=zform, data=demog,
    spot.check = TRUE, 
    ncheck = 10)



  ################
  # Check against truth
  ###############
  if (Sys.info()[["user"]] %in% "jmuschel") {
    x = load(file.path(outdir,
        paste0(adder, "Top_", runmod, "_Pvalues_df.Rda")))

  } else {
    x = load(file.path(rootdir, "CT_Pipeline",
      paste0(adder, "Top_", runmod, "_Pvalues_df.Rda")))

  }
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
  reg.smod = summary(reg.mod)

  zform2 = as.character(zform)
  zform2 = c(paste0("Y", zform2[1]), zform2[2])
  zform2 = as.formula(zform2)
  nox.mod = lm(formula=zform2, data=demog)
  nox.smod = summary(nox.mod)  

  # runmeas = c("r.squared", "adj.r.squared")
  runmeas = c("adj.r.squared")
  measure = "adj.r.squared"
  for (measure in runmeas){
    truth = reg.smod[[measure]]
    nox.truth = nox.smod[[measure]]
    perm = permute.mod[[measure]]

    # fit.gamma = function(x){
    #   mx = mean(x)
    #   n = length(x)
    #   ss = sum( (x - mx)^2 ) 
    #   ahat = mx^2 / ss
    #   thetahat = ss / (n * mx)
    #   return(list(shape=ahat, scale= thetahat))
    # }
    # gperm = fit.gamma(perm)
    fit.beta = function(x){
      m = mean(x)
      m2 = mean(x^2)
      u = m * (m - m2) / (m2 - m^2)
      v = (1 - m) * (m - m2) / (m2 - m^2)
      return(list(shape1=u, shape2= v))
    }    
    gbeta = fit.beta(perm)

    p = permute.mod$p
    ### p includes intercept
    n = permute.mod$n
    df = permute.mod$df

    rb = rbeta(10000, shape1= (p-1)/2, shape2= (n-p)/2)
    arb = 1 - (1-rb) * (n - 1) / (n - p)

    xlim = c(min(truth, perm), max(truth, perm))
    hist(perm, xlim=xlim, breaks=30, xlab=
      paste0("Permutation ", measure), 
      main=paste0("Permutation Distribution for ", measure, " on ",
      outcome, " using ", runmod, " HPR"))
    abline(v= truth)
    abline(v=nox.truth, col="red")
  }


  
  outfile = file.path(outdir, 
                      paste0(adder, "Permutation_Distribution.Rda"))  
  save(permute.mod, file=outfile)
  
}
dev.off()



# real.mod = fast_lm(orig.Y, X= as.matrix(perc), Z = zform, 
#   data=demog, spot.check=TRUE)
# real.r2 = real.mod$r.squared
# check.r2 = summary(reg.mod)$r.squared
# r.r2 = range(permute.mod$r.squared)
# xlim = c(min(check.r2, r.r2), max(check.r2, r.r2))
# hist(permute.mod$r.squared, xlim=xlim)
# abline(v= check.r2)
