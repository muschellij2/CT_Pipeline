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

fname = file.path(datadir, "Checked_Image_Filenames.csv")
check.imgs = read.csv(fname, as.is=TRUE)

fname = file.path(datadir, "SS_Check.txt")
fail = readLines(fname)
fail = strsplit(fail, " ")
fail = lapply(fail, function(x) x[!x %in% "-"])
fail = t(sapply(fail, function(x) x[1:2]))
fail = data.frame(fail, stringsAsFactors = FALSE)
colnames(fail) = c("img", "reason")
fail$img = gsub("[.]png", "", fail$img)
fail$stub = gsub("(.*)_SS.*", "\\1", fail$img)
fail$level = 3
fail$level[fail$reason %in% "eh"] = 1
fail$level[fail$reason %in% "CTA"] = 1

######################################
# Read in Craniotomy data
######################################
fname = file.path(datadir, "Craniotomy_Patients.csv")
crani = read.csv(fname, as.is=TRUE)
crani$id = paste0(floor(crani$patientName/1000), "-", 
	crani$patientName%% 1000)
crani$Time[ crani$Time == "" ] = "00:00"
crani$crani_date = paste(crani$Date_Time_ICH_Evacuation,
	crani$Time)
crani$crani_date = as.POSIXct(crani$crani_date, 
	format = "%d%b%Y %H:%M")

crani = crani[, c("id", "crani_date")]

######################################
# Get manual check data, 0 - bad, 0.25 - bad gantry tilt correct
# 0.5 - ok gantry tilt correct, but brain deleted
# 1 - good
######################################
check.imgs = check.imgs[ check.imgs$Good > 0, ]
check.imgs$stub = gsub("[.]png$", "", nii.stub(check.imgs$img))

######################################
# Convert date to POSIXct
######################################
ctdate_to_date = function(x){
  days = gsub("(.*)_(.*)", "\\1", x)
  year = substr(x=days, start=1, stop=4)
  mon = substr(x=days, start=5, stop=6)
  day = substr(x=days, start=7, stop=8)  
  
  times = gsub("(.*)_(.*)", "\\2", x)
  times[times == "NA"] = "0000"
  hour = substr(times, 1, 2)
  min = substr(times, 3, 4)
  add.day = rep(0, length=length(x))

  add.hour = (min == "60")
  min[add.hour] = "00"
  hour[add.hour] = sprintf("%02.0f", 
  	as.numeric(hour[add.hour]) + 1)

  #########################
  # Paste together and make a Date Object
  #########################
  time = paste0(hour, ":", min)
  dt = as.POSIXct(paste(days, time), format = "%Y%m%d %H:%M")
  dt = dt + add.day
  dt
}


rerun = FALSE

ids = list.files(basedir, pattern = "^\\d.*\\d$")
iddirs = file.path(basedir, ids)
ssdirs = file.path(iddirs, "Skull_Stripped")

iid <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(iid)) iid = 5

######################################
# read in volume data
######################################
# adder = ""
adder = "_nopresmooth"
bad.ids = NULL
all.df = NULL
for (iid in seq_along(ids)){
	id = ids[iid]
	iddir = iddirs[iid]
	ssdir = ssdirs[iid]
	outdir = file.path(iddir, "results")
	if (!file.exists(outdir)){
		dir.create(outdir)
	}
	outfile = file.path(outdir, paste0(
		"Skull_Strip_Volumes", adder, ".Rda"))
	if (file.exists(outfile)){
		load(outfile)
		all.df = rbind(all.df, df)
	} else {
		bad.ids = c(bad.ids, id)
	}
	# print(iid)
}


######################################
# subset data
######################################
all.df$stub = nii.stub(all.df$fname, bn=TRUE)
all.df$int = gsub(".*SS_(.*)_Mask.*", "\\1", all.df$stub)
all.df$stub = gsub("(.*)_SS_(.*)_Mask.*", "\\1", all.df$stub)

all.df = all.df[ all.df$stub %in% check.imgs$stub, ]

all.df$date = gsub(".*(\\d{8}_(\\d{4}|NA)).*", "\\1", all.df$stub)
all.df$date = gsub(".*(\\d{8}_(\\d{4}|NA)).*", "\\1", all.df$stub)

all.df$date = ctdate_to_date(all.df$date)

all.df$id = gsub("(.*)_\\d{8}.*", "\\1", all.df$stub)
all.df$id = gsub("(.*)__.*", "\\1", all.df$id)
all.df$id = gsub("(.*)__.*", "\\1", all.df$id)

all.df = all.df[ order(all.df$id, all.df$date), ]

all.df = ddply(all.df, .(id), function(x){
	x$basedate = x$date[1]
	x
})

all.df$from_base = difftime(all.df$date, 
	all.df$basedate, units="days")

all.df = merge(all.df, crani, by="id", all.x=TRUE, sort=FALSE)

all.df = merge(all.df, fail, by="stub", all.x=TRUE, sort=FALSE)

all.df$level[is.na(all.df$level)] = 0

all.df = all.df[ order(all.df$id, all.df$date), ]

all.df$before_crani = TRUE
all.df$before_crani[!is.na(all.df$crani_date)] = 
	all.df$date[!is.na(all.df$crani_date)]  <=
	all.df$crani_date[!is.na(all.df$crani_date)] 

all.df$site_number = substr(all.df$id, 1, 3)

all.df$img = file.path(dirname(dirname(all.df$fname)), 
	gsub("(.*)_SS_.*", "\\1.nii.gz", basename(all.df$fname)))

stopifnot(all(file.exists(all.df$img)))

outfile = file.path(resdir, 
	paste0("Longitudinal_Skull_Strip_Data", adder, ".Rda"))
save(all.df, file= outfile)

set.seed(20141103)
all.df$alpha = 0.2
ids = sample(all.df$id, size = 5)
all.df$alpha[ all.df$id %in% ids] = 0.8
######################################
# subset 0.01
######################################	
df = all.df[ grep("SS_0.01", all.df$fname, fixed=TRUE), ]

check.ids = c("120-376", "131-316", "133-417", "134-327")
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
gtrue = g + aes(y= truevol) + 
	ylab("Estimated Intracranial Volume")
gstrue = gtrue  + 
	geom_smooth(aes(group=1), se=FALSE, col="black", size=1)
pngname = file.path(resdir, "Intraclass_Correlation.png")
png(pngname, type="cairo", 
	res=600, height=7, width=7, units="in")
	gstrue
dev.off()

mod = lmer(truevol ~ 1 + (1 | id), data=df)
smod = summary(mod)

vc = VarCorr(mod)
tau.sq <- vc$id[1]
sigma.sq <- smod$sigma^2
tau.sq/(tau.sq+sigma.sq)

nocrani = df[df$before_crani,]

pngname = file.path(resdir, "Intraclass_Correlation_no_crani.png")
png(pngname, type="cairo", 
	res=600, height=7, width=7, units="in")
	gstrue %+% nocrani
dev.off()
mod = lmer(truevol ~ 1 + (1 | id), data=nocrani)
smod = summary(mod)

vc = VarCorr(mod)
tau.sq <- vc$id[1]
sigma.sq <- smod$sigma^2
tau.sq/(tau.sq+sigma.sq)

ddf = df[df$before_crani & 
	df$level == 0,]

d=data.frame(x1=c(0), x2=c(10), 
	y1=c(min(ddf$truevol)), 
	y2=max(ddf$truevol), truevol=mean(ddf$truevol),
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
	"Intraclass_Correlation_no_crani_check.png")
png(pngname, type="cairo", 
	res=600, height=7, width=7, units="in")
	gstrue %+% ddf
dev.off()

pngname = file.path(resdir, 
	"Intraclass_Correlation_no_crani_check_fill.png")
png(pngname, type="cairo")
	gs2 %+% ddf
dev.off()
mod = lmer(truevol ~ 1 + (1 | id), data=ddf)
smod = summary(mod)

vc = VarCorr(mod)
tau.sq <- vc$id[1]
sigma.sq <- smod$sigma^2
tau.sq/(tau.sq+sigma.sq)

save(ddf, mod, smod, vc, 
	tau.sq, sigma.sq, 
	file = file.path(resdir, "ICC_data.Rda"))


ddf10 = df[df$before_crani & df$level == 0 & df$from_base < 10,]



pngname = file.path(resdir, 
	"Intraclass_Correlation_no_crani_check_day10.png")
png(pngname, type="cairo", 
	res=600, height=7, width=7, units="in")
	gstrue %+% ddf10
dev.off()

mod = lmer(truevol ~ 1 + (1 | id), data=ddf10)
smod = summary(mod)

vc = VarCorr(mod)
tau.sq <- vc$id[1]
sigma.sq <- smod$sigma^2
tau.sq/(tau.sq+sigma.sq)

pngname = file.path(resdir, 
	"Intraclass_Correlation_no_crani_check_day10_black.png")
png(pngname, type="cairo", 
	res=600, height=7, width=7, units="in")
	g10 = (gbase + aes(group=id, alpha=alpha)) %+% ddf10
	g10 = g10 + aes(y= truevol) + 
	ylab("Estimated Intracranial Volume") + 
	guides(alpha=FALSE)
	print(g10)
dev.off()

# "173-312" is bad one
g = ggplot(df[ df$site_number %in% c("173", "232"), ], 
	aes(x= as.numeric(from_base), 
	y= truevol, colour=id)) + 
	geom_line() + facet_wrap(~ site_number)

# gz = g + aes(y= zvol)
# gz
# gthick = g + aes(y = thickvol)
# gthick
# g = ggplot(df, aes(x= as.numeric(from_base), 
# 	y= truevol, colour=id)) + 
# 	geom_line() +
# 	guides(colour=FALSE)	