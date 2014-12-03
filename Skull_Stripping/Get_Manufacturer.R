###################################################################
## This code is for longitudinal skull stripping volume estimate
##
## Author: John Muschelli
###################################################################
###################################################################
rm(list=ls())
library(fslr)
library(cttools)
library(plyr)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
rerun = FALSE

progdir = file.path(rootdir, "programs")
basedir = file.path(rootdir, "Registration")
datadir = file.path(basedir, "data")
resdir = file.path(basedir, "results")
outdir = file.path(resdir, "SS_Overlays")
if(!file.exists(outdir)){
	dir.create(outdir)
}
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")
rundir = file.path(progdir, "Skull_Stripping")

datadir = file.path(basedir, "data")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")
rundir = file.path(progdir, "Skull_Stripping")

########################################
# Read in Checked data
########################################
cuts = c("0.01", "0.1", "0.35")
checks = lapply(cuts, function(cut){
	fname = file.path(datadir, 
		paste0("Check_", cut, ".csv"))
	x = read.csv(fname, as.is=TRUE)	
	x$int = cut
	x
})
checks = do.call("rbind", checks)
checks$name = gsub("/", ":", checks$name)
checks$stub = gsub("(.*)_SS_0.*", "\\1", checks$name)
checks = checks[ order(checks$stub, checks$int), ]
checks = ddply(checks, .(stub), function(x){
	x$Gantry = max(x$Gantry)
	x$Crani = max(x$Crani)
	x
})

check.imgs = checks
rm(list="checks")

outfile = file.path(resdir, "Longitudinal_Skull_Strip_Data.Rda")
load(file= outfile)

smooth = FALSE
if (!smooth){
	all.df$fname = gsub("_Mask", "_nopresmooth_Mask", all.df$fname)	
}

all.df = all.df[ file.exists(all.df$hdr), ]
all.df$manu = NA
all.df = all.df[, c("id", "stub", "hdr", "manu")]
all.df = unique(all.df)
for (idf in seq(nrow(all.df))){
	load(all.df$hdr[idf])
	manu = unique(dcmtables[, "0008-0070-Manufacturer"])
	stopifnot(length(manu) <= 1)
	all.df$manu[idf] = manu
	print(idf)
}

hdrs = all.df
outfile = file.path(resdir, "Manufacturer_Data.Rda")
save(hdrs, file= outfile)