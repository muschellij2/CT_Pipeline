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
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
# source(file.path(progdir, "convert_DICOM.R"))
# source(file.path(progdir, "fslhd.R"))
basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

outdir = file.path(basedir, "results")

whichdir = "reoriented"

get.id = function(x){
	ss = strsplit(x, "_")
	ss = sapply(ss, head, 1)
	ss = gsub(".*(\\d\\d\\d-.*)", "\\1", ss)
	ss
}

id_to_pname = function(x){
	as.numeric(gsub("-", "", x))
}

nkeeps = c(1000, 2000, 3000, 0.05, 0.01, .001)
cn = c("With_Perc", "With_Clot", "Null", 
	"Nosex_With_Perc", "Nosex_With_Clot", "Nosex_Null", "nkeep")


demog = read.csv(file=file.path(basedir, "Demog_NIHSS_Mask.csv"), 
	stringsAsFactors=FALSE)
demog$Base_ICH_10 = demog$Diagnostic_ICH /10

measures  = c("adj.r.squared", "r.squared")
reses = vector("list", length=length(measures))
names(reses) = measures
aics = reses 

for (meas in measures){

	res = matrix(NA, ncol = length(cn), nrow=length(nkeeps))
	colnames(res) = cn
	rownames(res) = nkeeps

	aic = res

	ikeep = 1

	for (ikeep in seq_along(nkeeps)){
		
		nkeep = nkeeps[ikeep]
		# nkeep = .001
	# function(nkeep){
		outfile = file.path(outdir, 
			paste0("Top_", nkeep, "_Pvalues_df.Rda"))
		load(file=outfile)


		perc = colMeans(submat)
		check.ids = id_to_pname(get.id(names(perc)))
		stopifnot(all(demog$patientName == check.ids))

		demog$perc_ROI = perc * 100

		mod = lm(Enrollment_NIHSS_Total ~ perc_ROI + Age + Sex + 
			Base_ICH_10, data=demog)
		smod = summary(mod)
		nmod = lm(Enrollment_NIHSS_Total ~ Age + Sex + 
			Base_ICH_10, data=demog)
		snmod = summary(nmod)
		cmod = lm(Enrollment_NIHSS_Total ~ Age + Sex + 
			Base_ICH_10 + Lobar_BG_vs, data=demog)
		scmod = summary(cmod)


		anova(mod, nmod)
		res[ikeep, "With_Perc"] = smod[[meas]]
		res[ikeep, "With_Clot"] = scmod[[meas]]
		res[ikeep, "Null"] = snmod[[meas]] 

		aic[ikeep, "With_Perc"] = AIC(mod)
		aic[ikeep, "With_Clot"] = AIC(cmod)
		aic[ikeep, "Null"] = AIC(nmod)
		#################################
		# Run model without sex
		#################################

		mod = lm(Enrollment_NIHSS_Total ~ perc_ROI + Age  + 
			Base_ICH_10, data=demog)
		smod = summary(mod)
		nmod = lm(Enrollment_NIHSS_Total ~ Age + 
			Base_ICH_10, data=demog)
		snmod = summary(nmod)
		cmod = lm(Enrollment_NIHSS_Total ~ Age + 
			Base_ICH_10 + Lobar_BG_vs, data=demog)
		scmod = summary(cmod)


		anova(mod, nmod)
		res[ikeep, "Nosex_With_Perc"] = smod[[meas]]
		res[ikeep, "Nosex_With_Clot"] = scmod[[meas]]
		res[ikeep, "Nosex_Null"] = snmod[[meas]] 

		aic[ikeep, "Nosex_With_Perc"] = AIC(mod)
		aic[ikeep, "Nosex_With_Clot"] = AIC(cmod)
		aic[ikeep, "Nosex_Null"] = AIC(nmod)
		
		aic[ikeep, "nkeep"] = nkeep
		res[ikeep, "nkeep"] = nkeep
	}
	reses[[meas]] = res
	aics[[meas]] = aic
}


save(reses, measures, nkeeps, aics, 
	file=file.path(outdir, "Regress_ROI_NIHSS_Results.Rda"))
