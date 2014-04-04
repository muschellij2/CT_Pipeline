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
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
source(file.path(progdir, "convert_DICOM.R"))
source(file.path(progdir, "fslhd.R"))
basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

outdir = file.path(basedir, "results")

whichdir = "reoriented"

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
where = "NIHSS_Regress.R"

output = file.path(outdir, "Voxel_Info.Rda")
save(nvox, ncut, nuniq.rows, 
	dups, group, where, 
	file=output)

class(mat) = "numeric"
#### baseline NIHSSS ata 
nihss = read.csv(file.path(basedir, 
	"baseline_NIHSS.csv"), 
                 stringsAsFactors=FALSE)
nihss = nihss[ nihss$patientName %in% df$id, ]
nihss = nihss[ order(nihss$patientName), ]

demog = read.csv(file.path(basedir, 
	"All_180_FollowUp_wDemographics.csv"), 
                 stringsAsFactors=FALSE)
#### creating factors
demog$Sex = factor(demog$Gender, levels=c("Female", "Male"))
demog$Diagnostic_ICH = demog$ICH_Dx_10 *10
demog$Clot_Location_RC =  factor(demog$Clot_Location_RC, 
	levels=c("Putamen", "Lobar", "Globus Palidus", "Thalamus"))
stopifnot(all(!is.na(demog$Clot_Location_RC)))
demog$Lobar_BG_vs <- as.character(demog$Lobar_BG_vs)
demog$Lobar_BG_vs[demog$Lobar_BG_vs == "Not Lobar"] <- "Deep"
demog$Lobar_BG_vs <- factor(demog$Lobar_BG_vs, c("Lobar", "Deep"))

### keep subset
demog = demog[ demog$patientName %in% df$id, ]
dd = demog[ order(demog$patientName), ]
stopifnot(identical(nihss$nihss_total, dd$Enrollment_NIHSS_Total))

write.csv(demog, file=file.path(basedir, "Demog_NIHSS_Mask.csv"), 
	row.names=FALSE)


# mod = lm(Enrollment_NIHSS_Total ~ Clot_Location_RC, data=demog)
# mod = lm(Enrollment_NIHSS_Total ~ Clot_Location_RC + Age, 
# 	data=demog)

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
# mat = mat[1:1000, ]
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
	mod.1 = Enrollment_NIHSS_Total ~ X,
	mod.2 = Enrollment_NIHSS_Total ~ X + Age,
	mod.3 = Enrollment_NIHSS_Total ~ X + Sex,		
	mod.4 = Enrollment_NIHSS_Total ~ X + Diagnostic_ICH,		
	mod.5 = Enrollment_NIHSS_Total ~ X + Age + Sex + Diagnostic_ICH,
	mod.6 = Enrollment_NIHSS_Total ~ X
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



outfile = file.path(outdir, "Pvalue_Matrix.Rda")
save(res, file=outfile)