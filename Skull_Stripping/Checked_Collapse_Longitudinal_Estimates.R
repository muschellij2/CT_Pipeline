###################################################################
## This code is for longitudinal skull stripping volume estimate
##
## Author: John Muschelli
###################################################################
###################################################################
rm(list=ls())
library(plyr)
library(cttools)
library(fslr)
library(matrixStats)
library(ggplot2)
library(lme4)
library(ICC)
library(psychometric)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
basedir = file.path(rootdir, "Registration")
datadir = file.path(basedir, "data")
resdir = file.path(basedir, "results")
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
sm.checks = lapply(cuts, function(cut){
	fname = file.path(datadir, 
		paste0("Check_", cut, "_nopresmooth.csv"))
	x = read.csv(fname, as.is=TRUE)	
	x$int = cut
	x
})
checks = do.call("rbind", checks)
sm.checks = do.call("rbind", sm.checks)
checks = rbind(checks, sm.checks)
checks$name = gsub("/", ":", checks$name)
checks$smooth = !grepl("nopresmooth", checks$name)
checks$stub = gsub("(.*)_SS_0.*", "\\1", checks$name)
checks = checks[ order(checks$stub, -checks$smooth, checks$int), ]
checks = ddply(checks, .(stub), function(x){
	x$Gantry = max(x$Gantry)
	x$Crani = max(x$Crani)
	x
})

checks = ddply(checks, .(stub), function(x){
	x$anyNA = any(is.na(x$Good))
	x
})

check.imgs = checks
rm(list="checks")


#############################
check.imgs$name = gsub("[.]png$", "", check.imgs$name)

outfile = file.path(resdir, "Longitudinal_Skull_Strip_Data.Rda")
load(file= outfile)

########################################
# Merge Volume data
########################################
all.df$reason = all.df$level = NULL

all.df$name = nii.stub(all.df$fname, bn=TRUE)

all.df$all = 1
check.imgs$checked = 1

sm.checks = check.imgs[ !check.imgs$smooth, ]

check.imgs = check.imgs[ check.imgs$smooth, ]

all.df = merge(all.df, check.imgs, all = TRUE, 
	by = c("name", "int", "stub"))


########################################
# Double check all are checked
########################################
stopifnot(!any(is.na(all.df$all)))
stopifnot(!any(is.na(all.df$checked)))

all.df$all = all.df$checked = NULL
all.df$reason = all.df$level = NULL
xall.df = all.df

outfile = file.path(resdir, 
	"Longitudinal_Skull_Strip_Data_nopresmooth.Rda")
load(file= outfile)

all.df$all = 1
sm.checks$checked = 1

all.df$int = gsub("_nopresmooth", "", all.df$int)
all.df$name = nii.stub(all.df$fname, bn=TRUE)

all.df = merge(all.df, sm.checks, all.y = TRUE, 
	by = c("name", "int", "stub"))

###################################
# This needs to be removed next round
###################################
all.df = all.df[!is.na(all.df$all) & 
	!is.na(all.df$checked), ]

stopifnot(!any(is.na(all.df$all)))
stopifnot(!any(is.na(all.df$checked)))

all.df$all = all.df$checked = NULL
all.df$reason = all.df$level = NULL

xall.df = xall.df[ xall.df$img %in% all.df$img, ]
xall.df = xall.df[, colnames(all.df)]

all.df = rbind(xall.df, all.df)

all.df = all.df[!all.df$anyNA, ]
all.df$anyNA = NULL

test = all.df[ all.df$int == "0.01" &
	all.df$smooth == TRUE, ]

###########################
# Getting N's for paper
###########################
outfile = file.path(resdir, "Manufacturer_Data.Rda")
load(file= outfile)

test = merge(test, hdrs[, c("hdr", "manu")], 
	all.x=TRUE, all.y=FALSE)
stopifnot(all(!is.na(test$manu)))

all.manu.tab = table(test$manu)


total.N = nrow(test)
N.crani = sum(test$Crani >= 1)
N.gantry = sum(test$Gantry >= 1 & test$Crani < 1)

test = test[ test$Gantry < 1 & test$Crani < 1, ]

nppt = ddply(test, .(id), nrow)
n.per.pt = quantile(nppt$V1, probs = c(0, 0.25, 0.5, 0.75, 1))
n.per.pt = c(n.per.pt, mean=mean(nppt$V1), sd=sd(nppt$V1))

n.ctr.icc = length(unique(test$site_number))
manu.tab = table(test$manu)
total.Npt = length(unique(test$id))


all.df = merge(all.df, hdrs[, c("hdr", "manu")], 
	all.x=TRUE, all.y=FALSE)
stopifnot(all(!is.na(all.df$manu)))


n.crani = sum(all.df$Crani < 1)
nocrani = all.df[all.df$Crani < 1, ]



check.crani = nocrani[ ! nocrani$before_crani ,]
check.crani = unique(
	check.crani[, c("img", "date", "crani_date", "id")])
check.crani = ddply(check.crani, .(id), function(x){
	x = x[ which(x$date == max(x$date)),]
	x
})

########################
## these have "craniotomies", but aren't big or apparent -keep
########################
crani.ids = c("100-318", "102-349", "223-355", 
	"225-507", "265-398")
stopifnot(all(check.crani$id %in% crani.ids))


fail.rates = ddply(nocrani, .(int, smooth), function(x){
	c(Gantry= mean(x$Gantry ==1), 
		Bad.NoGant.pct = mean(x[ x$Gantry < 1, "Good"] < 1),
		Good.NoGant.pct = mean(x[ x$Gantry < 1, "Good"] == 1),
		Bad.NoGant.sum = sum(x[ x$Gantry < 1, "Good"] < 1),
		Good.NoGant.sum = sum(x[ x$Gantry < 1, "Good"] == 1),
		N.NoGant = sum(x$Gantry < 1),
		N = nrow(x)
		)
})



scan.fail.rates = ddply(nocrani, .(int, smooth, manu), function(x){
	c(Gantry= mean(x$Gantry ==1), 
		Bad.NoGant.pct = mean(x[ x$Gantry < 1, "Good"] < 1),
		Good.NoGant.pct = mean(x[ x$Gantry < 1, "Good"] == 1),
		Bad.NoGant.sum = sum(x[ x$Gantry < 1, "Good"] < 1),
		Good.NoGant.sum = sum(x[ x$Gantry < 1, "Good"] == 1),
		N.NoGant = sum(x$Gantry < 1),
		N = nrow(x)
		)
})

sfail_0.01 = scan.fail.rates[ scan.fail.rates$int == 0.01 & 
	scan.fail.rates$smooth== TRUE,]

###############################
# Sampling Ids
###############################
set.seed(20141119)
all.df$alpha = .4
ids = sample(all.df$id, size = 10)
all.df$alpha[ all.df$id %in% ids] = .8
# all.df$alpha = factor(all.df$alpha, levels= c("Low", "High"))

nocrani = all.df[all.df$Crani < 1, ]


######################################
# subset each int
######################################	

run.int = "0.01"
runvol = "truevol"
smooth = TRUE
runvols =  c("truevol", "thickvol")
ints = c("0.01", "0.35", "0.1")
smooths = c(TRUE, FALSE)
icc.vals = expand.grid(runvol = runvols,
	smooth = smooths,
	int = ints, stringsAsFactors=FALSE)
icc.vals$good.ICC.novar = icc.vals$good.ICC = icc.vals$ICC = NA
iscen = 1
icc.mods = vector(mode="list", length= nrow(icc.vals))
icc.smods = icc.mods
icc.cis = icc.mods

icc.allcis = icc.allsmods = icc.allmods = icc.mods


for (iscen in seq(nrow(icc.vals))) {

	runvol = icc.vals$runvol[iscen]
	run.int = icc.vals$int[iscen]
	smooth = icc.vals$smooth[iscen]

	df = nocrani[ nocrani$int == run.int 
		& nocrani$smooth == smooth, ]
	df$runvol = df[, runvol]

	addsmooth = ""
	if (!smooth){
		addsmooth = "_nopresmooth"
	}
	# check.ids = c("120-376", "131-316", "133-417", "134-327")
	tsize = 16

	gbase = ggplot(df, aes(x= as.numeric(from_base))) + 
		geom_line() +
		guides(colour=FALSE)  + 
		xlab("Days from Baseline Scan") +
		ggtitle("Estimate of Intracranial Volume Over Time") + 
	  theme(legend.position = c(.5, .5),
	        legend.background = element_rect(fill="transparent"),
	        legend.key = element_rect(fill="transparent", 
	                                  color="transparent"),
	        legend.text = element_text(size=tsize+2), 
	        legend.title = element_text(size=tsize),
	        title = element_text(size=tsize),
	        strip.text = element_text(size = tsize+4),
	        axis.text  = element_text(size=tsize-2))
	g = gbase + aes(colour=id)
	gtrue = g + aes(y= runvol) + 
		ylab("Estimated Intracranial Volume (cc)")
	gstrue = gtrue  + 
		geom_smooth(aes(group=1), se=FALSE, col="black", size=1)
	pngname = file.path(resdir, 
		paste0("Intraclass_Correlation_", run.int, 
			"_", runvol, addsmooth, ".png"))
	png(pngname, type="cairo", 
		res=600, height=7, width=7, units="in")
		gstrue
	dev.off()

	mod = lmer(runvol ~ 1 + (1 | id), data=df)
	smod = summary(mod)

	vc = VarCorr(mod)
	tau.sq <- vc$id[1]
	sigma.sq <- smod$sigma^2

	icc.vals$ICC[iscen] = tau.sq/(tau.sq+sigma.sq)
	icc.allmods[[iscen]] = mod	
	icc.allsmods[[iscen]] = smod
	icc.allcis[[iscen]] = ICCest("id", "runvol", data=df)

	####################################
	# Subset only good data
	####################################

	ddf = df[ df$Good == 1,]

	d=data.frame(x1=c(0), x2=c(10), 
		y1=c(min(ddf$runvol)), 
		y2=max(ddf$runvol), runvol=mean(ddf$runvol),
		from_base=0)

	grect = geom_rect(data=d, 
		mapping=aes(NULL, NULL, 
			xmin=x1, xmax=x2, ymin=y1, ymax=y2), 
		color="black", alpha=0.5, fill = "blue") 
	gs2 = gstrue + geom_rect(data=d,
		aes(xmin=x1, xmax=x2, ymin=-Inf, ymax=Inf,
			data=d),
		color=NA, alpha=0.5, fill = "blue") 


	pngname = file.path(resdir, 
		paste0("Intraclass_Correlation_no_crani_check_", run.int, 
			"_", runvol, addsmooth, ".png"))
	png(pngname, type="cairo", 
		res=600, height=7, width=7, units="in")
		gstrue %+% ddf
	dev.off()

	pngname = file.path(resdir, 
		paste0("Intraclass_Correlation_no_crani_check_fill_", 
			run.int, 
			"_", runvol, addsmooth, ".png"))
	png(pngname, type="cairo")
		gs2 %+% ddf
	dev.off()
	mod = lmer(runvol ~ 1 + (1 | id), data=ddf)
	smod = summary(mod)

	vc = VarCorr(mod)
	tau.sq <- vc$id[1]
	sigma.sq <- smod$sigma^2

	icc.vals$good.ICC[iscen] = tau.sq/(tau.sq+sigma.sq)
	icc.mods[[iscen]] = mod	
	icc.smods[[iscen]] = smod
	icc.cis[[iscen]] = ICCest("id", "runvol", data=ddf)

	####################################
	# Subset No variable slice thickness scans
	####################################

	ddf.novar = ddf[ !ddf$varslice,]

	pngname = file.path(resdir, 
		paste0("Intraclass_Correlation_no_crani_check_novarslice_", 
			run.int, 
			"_", runvol, addsmooth, ".png"))
	png(pngname, type="cairo", 
		res=600, height=7, width=7, units="in")
		gstrue %+% ddf.novar
	dev.off()

	mod = lmer(runvol ~ 1 + (1 | id), data=ddf.novar)
	smod = summary(mod)

	vc = VarCorr(mod)
	tau.sq <- vc$id[1]
	sigma.sq <- smod$sigma^2

	icc.vals$good.ICC.novar[iscen] = tau.sq/(tau.sq+sigma.sq)



	####################################
	# Saving Results
	####################################
	save(ddf, mod, smod, vc, 
		tau.sq, sigma.sq, 
		file = file.path(resdir, 
		paste0("ICC_data_", 
			run.int, 
			"_", runvol, addsmooth, ".Rda"))
		)


	####################################
	# Keep only Acute phase < 10 days post baseline
	####################################
	ddf10 = ddf[ ddf$from_base < 10,]



	pngname = file.path(resdir, 
		paste0("Intraclass_Correlation_no_crani_check_day10_", 
			run.int, 
			"_", runvol, addsmooth, ".png"))
	png(pngname, type="cairo", 
		res=600, height=7, width=7, units="in")
		gstrue %+% ddf10
	dev.off()

	mod = lmer(runvol ~ 1 + (1 | id), data=ddf10)
	smod = summary(mod)

	vc = VarCorr(mod)
	tau.sq <- vc$id[1]
	sigma.sq <- smod$sigma^2
	tau.sq/(tau.sq+sigma.sq)

	pngname = file.path(resdir, 
		paste0(
			"Intraclass_Correlation_no_crani_check_day10_black_", 
			run.int, 
			"_", runvol, addsmooth, ".png"))
	png(pngname, type="cairo", 
		res=600, height=7, width=7, units="in")
		g10 = (gbase + aes(group=id, alpha=alpha)) %+% ddf10
		g10 = g10 + aes(y= runvol) + 
		ylab("Estimated Intracranial Volume (cc)") + 
		guides(alpha=FALSE) + geom_smooth(aes(group=1), 
			se=FALSE, col="blue", size=1) 
		print(g10)
	dev.off()

	print(iscen)
	# "173-312" is bad one
	# 134-412
	# sites = c("173", "232")
} # iscen

good.cis = cbind(
	good.lower=sapply(icc.cis, `[[`, "LowerCI"),
	good.est=sapply(icc.cis, `[[`, "ICC"),
	good.upper=sapply(icc.cis, `[[`, "UpperCI"))

cis = cbind(
	all.lower=sapply(icc.allcis, `[[`, "LowerCI"),
	all.est=sapply(icc.allcis, `[[`, "ICC"),
	all.upper=sapply(icc.allcis, `[[`, "UpperCI"))

icc.vals = cbind(icc.vals, cis, good.cis)

save(icc.vals, icc.mods, icc.cis, icc.smods,
	total.N, N.crani, N.gantry, total.Npt,
	all.manu.tab, manu.tab,
	n.ctr.icc,
	n.per.pt,
	fail.rates,
	file = file.path(resdir, 
	paste0("ICC_Results.Rda"))
	)

sites = "134"
g = ggplot(ddf[ ddf$site_number %in% sites, ], 
	aes(x= as.numeric(from_base), 
	y= runvol, colour=id)) + 
	geom_line() + facet_wrap(~ site_number)


g = ggplot(ddf[ ddf$id %in% "134-412", ], 
	aes(x= as.numeric(from_base), 
	y= runvol)) + 
	geom_line()

# gz = g + aes(y= zvol)
# gz
# gthick = g + aes(y = thickvol)
# gthick
# g = ggplot(df, aes(x= as.numeric(from_base), 
# 	y= truevol, colour=id)) + 
# 	geom_line() +
# 	guides(colour=FALSE)	