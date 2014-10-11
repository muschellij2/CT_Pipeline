rm(list=ls())
library(cttools)
library(plyr)
library(ggplot2)
library(reshape2)
homedir = "/Applications"
rootdir = "~/CT_Registration"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident/Registration"
}
basedir = file.path(rootdir, "data")
outdir = file.path(rootdir, "results")
load(file.path(outdir, "Voxel_Matrix.Rda"))

ncut = 10
N = ncol(mat)
mat = mat[rs > ncut & rs < (N - 10), ]
probs = rowMeans(mat)

X = t(mat)
n.sig = 20
set.seed(20141006)
sig = sample(nrow(mat), size = n.sig)

x = X[,sig ]
true.pct = rowMeans(x)
betas = seq(0, 2, by=0.25)
ibeta <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(ibeta)) ibeta = 2
beta = betas[ibeta]

B = 1000
E = matrix(rnorm(N*B), nrow=N, ncol = B)
BX = matrix(beta* true.pct, nrow=N, ncol = B)
Y = BX + E

# true.mod = lm(Y ~ true.pct)
# true.smod = summary(true.mod)
# true.r2 = true.smod$r.squared

# Y =  rnorm(N)


###
# Statistic - increase in r^2
# If data has no bearing on Y, then if you get the ROI pct 
# (many times), you get distribution fo R^2.
# What do we compare that to?
###

mods = fast_lm(Y = Y, X = X, ncheck = 20)

cuts = c(0.001, 0.01, 0.05, 1000, 2000, 3000)
all.r2 = matrix(NA, nrow=length(cuts), ncol= B)
cut.corr = N.thresh = all.r2
pvals = mods$p.val

corr.sig = vector(length=length(cuts), mode= "list")
corrs = matrix(NA, ncol= B, nrow=length(cuts))
iicut = 1


for (iicut in seq_along(cuts)){
	#### grab the ROIs
	icut = cuts[iicut]
	if (icut < 1){
		ind = alply(pvals, .margins = 2, function(x) {
			which(x <= icut)
		}, .progress = "text")
	} else {
		ind = alply(pvals, .margins = 2, function(x) {
			r = rank(x)
			ind = which(r <= icut)			
		}, .progress = "text")
	}
	######################## 
	# If no percentages, just use 1 - will become unstable
	########################
	threshes = laply(ind, length)
	N.thresh[iicut, ] = threshes

	corr.sig[[iicut]] = laply(ind, function(x){
		sig %in% x
	})

	cut.corr[iicut,] = laply(ind, function(x){
		m = mean(x %in% sig)
		m[is.nan(m)] = 0
		m
	})

	pcts = t(laply(ind, function(x){
		if (length(x) == 0){
			return(rep(0, N))
		}
		roi = X[, x, drop=FALSE]
		pct = rowMeans(roi)
	}, .progress = "text"))
	corrs[iicut,] = cor(pcts, true.pct)	
	pb = txtProgressBar(max=B, style=3)
	for (iB in seq(B)){
		mod = lm(Y[,iB] ~ pcts[,iB])
		smod = summary(mod)
		all.r2[iicut, iB] = smod$r.squared
		setTxtProgressBar(pb, iB)
	}
	close(pb)
	print(icut)
}

rownames(all.r2) = cuts
r2.df = melt(all.r2)
colnames(r2.df) = c("Cutoff", "iteration", "r2")
r2.df$Cutoff = round(r2.df$Cutoff , 5)

sig.df = t(sapply(corr.sig, rowMeans))
rownames(sig.df) = cuts
sig.df = melt(sig.df)
colnames(sig.df) = c("Cutoff", "iteration", "Corr_Include")
r2.df$Corr_Include = sig.df$Corr_Include

rownames(corrs) = cuts
cor.df = melt(corrs)
colnames(cor.df) = c("Cutoff", "iteration", "corr")
r2.df$corr = cor.df$corr

rownames(N.thresh) = cuts
N.df = melt(N.thresh)
colnames(N.df) = c("Cutoff", "iteration", "n_thresh")
r2.df$n_thresh = N.df$n_thresh

rownames(cut.corr) = cuts
ccut.df = melt(cut.corr)
colnames(ccut.df) = c("Cutoff", "iteration", "Included_Correct")
r2.df$Included_Correct = ccut.df$Included_Correct

save(all.r2, r2.df, beta, N.thresh, pvals, corrs, corr.sig,
	file = file.path(outdir, 
	paste0("Simulation_Results_Beta_", beta, ".Rda")))

# pdfname = file.path(outdir, 
# 	paste0("Simulation_Results_Beta_", beta, ".pdf"))
# pdf(pdfname)
# 	g = ggplot(r2.df, aes(x=r2)) + geom_histogram() + 
# 	facet_wrap(~ Cutoff)
# 	print(g)
# dev.off()


