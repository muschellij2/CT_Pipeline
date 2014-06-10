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

vox.nkeeps = rep(NA, length(nkeeps))
mods = list()
for (meas in measures){

	res = matrix(NA, ncol = length(cn), nrow=length(nkeeps))
	colnames(res) = cn
	# rownames(res) = nkeeps

	epic = aic = res

	ikeep = 1

	pdfname = file.path(outdir, 
		paste0("Regress_ROI_", outcome, "_Plots.pdf"))
	pdf(pdfname)
	for (ikeep in seq_along(nkeeps)){
		
		nkeep = nkeeps[ikeep]
		# nkeep = .001
	# function(nkeep){
		outfile = file.path(outdir, 
			paste0(adder, "Top_", nkeep, "_Pvalues_df.Rda"))
		load(file=outfile)
		vox.nkeeps[ikeep] = nkeep

		perc = colMeans(submat)
		check.ids = id_to_pname(get.id(names(perc)))
		stopifnot(all(demog$patientName == check.ids))

		demog$perc_ROI = perc * 100
		### divide by 10 for 10% increases
		demog$perc_ROI  = demog$perc_ROI / 10 

		g = ggplot(aes(y=Y, x=perc_ROI), data=demog) + 
			ggtitle(paste0("Number of Voxels: ", nkeep)) +
			geom_point() + geom_smooth() + 
			geom_smooth(se=FALSE, method="lm", col="red") +
			xlab("Percent HPR Engagement") + 
			ylab(outcome)
		print(g)

		if ( ((outcome == "GCS" & nkeep == 1000) |
			(outcome == "NIHSS" & nkeep == 19047))
			& meas == measures[1]
			){
			pngname = file.path(outdir, 
				paste0("Regress_ROI_", outcome, "_Best_Model.png"))
			png(pngname, res = 600, height=7, width=7, units = "in")
				g = ggplot(aes(y=Y, x=perc_ROI), data=demog) + 
					ggtitle(paste0(outcome, " Score-HPR Coverage Relationship")) +
					geom_point() + geom_smooth(se=FALSE,
						aes(colour="LOESS"), size=1.5) + 
					geom_smooth(se=FALSE, method="lm", 
						aes(colour="Linear Model"), size=1.5) +
					xlab(paste0("Percent HPR Engagement ", 
						"with Best Model Fit (V=", nkeep, ")")) + 
					ylab(paste0(outcome, " Score"))
				g = g + scale_colour_manual("", 
					values = c("LOESS" = "blue", "Linear Model" = "red"))
				g = g + guides(colour = guide_legend(override.aes = list(size=2)))
				g = g + theme(axis.title = element_text(size = rel(1.3)),
					axis.text = element_text(size=rel(1.1)),
					plot.title = element_text(size=rel(1.5))
					)

				g = g+ theme(legend.position=c(.5, .075), 
					legend.background = element_rect(fill = "transparent"),
					legend.key.size = unit(1, "cm"),
					legend.text = element_text(size = 20),
					legend.key = element_rect(fill = "transparent", 
						colour = "transparent"))
				print(g)
			dev.off()
		}

		mod = lm(Y ~ Age + Sex + 
			Base_ICH_10 + perc_ROI, data=demog)
		mods[[ikeep]] = mod
		smod = summary(mod)
		nmod = lm(Y ~ Age + Sex + 
			Base_ICH_10, data=demog)
		snmod = summary(nmod)
		keep.cmod = cmod = lm(Y ~ Age + Sex + 
			Base_ICH_10 + LOC, data=demog)
		scmod = summary(cmod)


		anova(mod, nmod)
		res[ikeep, "With_Perc"] = smod[[meas]]
		res[ikeep, "With_Clot"] = scmod[[meas]]
		res[ikeep, "Null"] = snmod[[meas]] 
		res[ikeep, "N_Gr0"] = sum(demog$perc_ROI > 0)

		aic[ikeep, "With_Perc"] = AIC(mod)
		aic[ikeep, "With_Clot"] = AIC(cmod)
		aic[ikeep, "Null"] = AIC(nmod)
		aic[ikeep, "N_Gr0"] = sum(demog$perc_ROI > 0)
  
	  	epic[ikeep, "With_Perc"] = AIC(mod, k=4)
	  	epic[ikeep, "With_Clot"] = AIC(cmod, k=4)
	  	epic[ikeep, "Null"] = AIC(nmod, k=4)  
		epic[ikeep, "N_Gr0"] = sum(demog$perc_ROI > 0)
		#################################
		# Run model without sex
		#################################

		mod = lm(Y ~ perc_ROI + Age  + 
			Base_ICH_10, data=demog)
		smod = summary(mod)
		nmod = lm(Y ~ Age + 
			Base_ICH_10, data=demog)
		snmod = summary(nmod)
		cmod = lm(Y ~ Age + 
			Base_ICH_10 + LOC, data=demog)
		scmod = summary(cmod)


		anova(mod, nmod)
		res[ikeep, "Nosex_With_Perc"] = smod[[meas]]
		res[ikeep, "Nosex_With_Clot"] = scmod[[meas]]
		res[ikeep, "Nosex_Null"] = snmod[[meas]] 

		aic[ikeep, "Nosex_With_Perc"] = AIC(mod)
		aic[ikeep, "Nosex_With_Clot"] = AIC(cmod)
		aic[ikeep, "Nosex_Null"] = AIC(nmod)
		
	  	epic[ikeep, "Nosex_With_Perc"] = AIC(mod, k=4)
	  	epic[ikeep, "Nosex_With_Clot"] = AIC(cmod, k=4)
	  	epic[ikeep, "Nosex_Null"] = AIC(nmod, k=4)  
  
		aic[ikeep, "nkeep"] = nkeep
	  	epic[ikeep, "nkeep"] = nkeep
	  	res[ikeep, "nkeep"] = nkeep

		aic[ikeep, "pval"]= pval
	  	epic[ikeep, "pval"]= pval
    	res[ikeep, "pval"] = pval		
	}
	dev.off()
	reses[[meas]] = res
	aics[[meas]] = aic
	epics[[meas]] = epic
}

aics = aics[[1]]
epics = epics[[1]]

loc.tab = sort(table(demog$LOC))
loc.ptab = round(prop.table(loc.tab)*100, 1)
loc.levs = names(loc.tab)


save(reses, measures, nkeeps, aics, epics, vox.nkeeps,
	loc.levs, loc.tab, loc.ptab, mods, keep.cmod,
	file=file.path(outdir, 
		paste0("Regress_ROI_", outcome, "_Results.Rda"))
	)
