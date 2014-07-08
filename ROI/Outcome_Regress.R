#################################
# Performs voxel-wise models for age/gender/volume
#################################
rm(list=ls())
library(oro.nifti)
library(plyr)
library(limma)
library(microbenchmark)
library(abind)
library(ggplot2)
library(data.table)
library(fslr)
library(cttools)
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
outcome = "NIHSS"
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

#### getting unique rows and duplications
dups = duplicated(mat)
# cn = paste0("V", 1:ncol(mat))
cn = colnames(mat)
mm = cbind(mat, id = 1:nrow(mat))
dt = data.table(mm, key=cn)
urows = unique(dt)
urows[, c('group') := 1:nrow(urows)]
urows[, c('id') := NULL]
xx = dt[urows]

group = xx[, list(id, group)]
setkey(group, 'id')

nuniq.rows = nrow(urows)
nvox = nrow(mat)
where = "Outcome_Regress.R"

output = file.path(outdir, paste0(adder, "Voxel_Info.Rda"))
save(nvox, ncut, nuniq.rows, 
	dups, group, where, all.nvox,
	file=output)

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


rm(list=c("nihss", "dd"))

cn= c("X", "Age", "SexMale")
### add 1 for global
ncn = length(cn) + 1
bnames = c("Xbeta")
nbeta = length(bnames)
nmod = 6
nlm = 5
# mat = mat[1:100, ]
res = array(NA, dim = c(nrow(mat), ncn+nbeta, nmod), 
	dimnames = list(
	seq(nrow(mat)), 
	c("global", cn, bnames), 
	paste("mod.", 1:nmod, sep='') )
	)

# function to get p-values from model
	t.pval = function(x)  coef(x)[, "Pr(>|t|)"]
	get.pval = function(mods){
		smods = llply(mods, summary)
		pval = laply(smods, function(x){
			fs = x$fstatistic
			f.pval = pf(fs[1L], df1=fs[2L], df2=fs[3L], 
				lower.tail=FALSE)
			tpval = t.pval(x)[cn]
			return(c(f.pval, tpval))
		}, .progress = "text")
		return(pval)
	}
## function to get beta values
	get.beta = function(mods){
		betas = laply(mods, function(x){
			x$coefficients['X']
		})
		return(betas)
	}

	ncn.seq = 1:ncn 
	nbeta.seq = (ncn+1):(ncn+nbeta)

### formulas
forms = list(
	mod.1 = Y ~ X,
	mod.2 = Y ~ X + Age,
	mod.3 = Y ~ X + Sex,		
	mod.4 = Y ~ X + Diagnostic_ICH,		
	mod.5 = Y ~ X + Age + Sex + Diagnostic_ICH,
	mod.6 = Y ~ X
	)
	

	imod = 1;
	demog$X = NA


	for (imod in 1:nlm){
		### make models
		mods = alply(mat, 1, .fun = function(X){
			demog$X = X
			mod = lm(formula=forms[[imod]], data=demog)
		}, .progress = "text")

		### save results
		print(paste0("Model ", imod, " Made"))
		res[, ncn.seq, imod] = get.pval(mods)
		res[, nbeta.seq, imod] = get.beta(mods)
		print(paste0("Results ", imod, " Made"))
	}

	if (nmod > nlm){
		for (imod in seq(from=nlm+1, to=nmod)) {
			### make models
			wt.pvals = aaply(mat, 1, .fun = function(X) {
				demog$X = X
				wt = wilcox.test(forms[[imod]], 
					data=demog, exact=FALSE)
				pval = wt$p.value
			}, .progress = "text")
			print(paste0("Wilcoxon ", imod, " Made"))

			### save results
			mat = matrix(NA, nrow=nrow(res), ncol= length(ncn.seq))
			mat[,2] = wt.pvals
			res[, ncn.seq, imod] = mat		
			print(paste0("Results ", imod, " Made"))
		}
	}



outfile = file.path(outdir, paste0(adder, "Pvalue_Matrix.Rda"))
save(res, file=outfile)