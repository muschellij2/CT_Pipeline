rm(list=ls())
library(cttools)
library(plyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)
homedir = "/Applications"
rootdir = "~/CT_Registration"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident/Registration"
}
basedir = file.path(rootdir, "data")
outdir = file.path(rootdir, "results")

betas = seq(0, 2, by=0.25)
cuts = c(0.001, 0.01, 0.05, 1000, 2000, 3000)


files = file.path(outdir, 
	paste0("Simulation_Results_Beta_", betas, ".Rda"))
all.df = NULL

for (ifile in seq_along(files)){
	file = files[ifile]
	load(file)
	r2.df$beta = beta
	all.df = rbind(all.df, r2.df)
	print(ifile)
	# all.r2, r2.df, beta, N.thresh, pvals, corrs, corr.sig,
}

head(all.df)


pdfname = file.path(outdir, 
	paste0("Simulation_Results.pdf"))
pdf(pdfname)
g = ggplot(all.df, aes(x=r2)) + geom_histogram() + 
	facet_wrap(Cutoff ~ beta)
print(g)


g = ggplot(all.df, aes(x=r2, colour = factor(beta))) + 
	geom_density() + 
	facet_wrap( ~ Cutoff)
print(g)

cols = c("red", rev(brewer.pal(9, "Blues")))
colfunc = seq_gradient_pal(low="purple",
	high = "deepskyblue")
cols = c("red", colfunc(1:9/9))
g = ggplot(all.df, aes(x=r2, colour = factor(beta))) + 
	geom_line(stat="density") + 
	facet_wrap( ~ Cutoff) + 
	scale_colour_manual(values = cols)
print(g + theme(legend.position = "bottom"))


cols = c("yellow",  "orange", "red",
	"green", "deepskyblue", "blue")
g = ggplot(all.df, aes(x=r2, colour = factor(Cutoff))) + 
	geom_line(stat="density") + 
	facet_wrap( ~ beta) + 
	scale_colour_manual(values = cols)
print(g + theme(legend.position = "bottom"))

dev.off()