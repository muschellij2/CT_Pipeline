#################################
# Performs voxel-wise models for age/gender/volume
#################################
rm(list=ls())
library(cttools)
library(matrixStats)
library(vows)
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

outdir = file.path(basedir, "results")

whichdir = "reoriented"
outcome = "GCS"
adder = paste0(outcome, "_")
if (outcome == "NIHSS"){
	adder = ""
}

#### load voxel data
outfile = file.path(outdir, "Voxel_Matrix.Rda")
load(file=outfile )

#### keeping if over 10 people have ICH in that locaiton
ncut = 10
all.nvox = sum(rs > 0)
mat = mat[rs > ncut, ]

class(mat) = "numeric"

demog = read.csv(file=file.path(basedir, "Demog_NIHSS_Mask.csv"), 
	stringsAsFactors=FALSE)

if (outcome == "GCS") {
	demog$Y = demog$Enrollment_GCS_Add
} else if (outcome == "NIHSS"){
	demog$Y = demog$Enrollment_NIHSS_Total
} else {
	stop(paste0("Outcome ", outcome, " not implemented"))
}


## table of ICES vs. MISTIE
nbreak = table(demog$patientName %% 1000 > 500)
sitetab = table(floor(demog$patientName/1000))
nsite = length(sitetab)


### formulas
forms = list(
	mod.1 = Y ~ X,
	mod.2 = Y ~ X + Age,
	mod.3 = Y ~ X + Sex,		
	mod.4 = Y ~ X + Diagnostic_ICH,		
	mod.5 = Y ~ X + Age + Sex + Diagnostic_ICH,
	mod.6 = Y ~ X
	)

rm(list=c("nihss", "dd"))
Y = demog$Y
Z = model.matrix(~ Age + Sex + Diagnostic_ICH, data= demog)
cc = complete.cases(Z) & complete.cases(Y)
Y = Y[cc]
Z = Z[cc,, drop=FALSE]
Zhat = Z %*% solve( t(Z) %*% Z ) %*% t(Z)
Yhat = Zhat %*% Y
R_Y = Y - Yhat

X = mat[, cc]
tX = t(X)
Xhat = Zhat %*% tX
R_tX  = tX - Xhat
rxvars = colVars(R_tX)
rxss = colSums(R_tX^2)
R_X = t(R_tX)

cov_RXY = R_X %*% R_Y
beta_x = cov_RXY / rxss