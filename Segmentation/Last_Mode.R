rm(list=ls())
library(fslr)
library(VGAM)
library(scales)
library(WhiteStripe)
x = readNIfTI("100-318_20070723_0957_CT_3_CT_Head-_Masked_SS.nii.gz")
roi = readNIfTI("100-318_20070723_0957_CT_3_CT_Head-ROI.nii.gz")

df = data.frame(val = c(x))
df$roi = c(roi)
sum(df$roi[ df$val < 20 | df$val > 100])
df = df[ df$val >= 20 & df$val <= 100, ]
s = sample(nrow(df), size=2e4)
samp = df[s,]


fit2 <- vglm(val ~ 1, mix2normal(eq.sd = TRUE, qmu=c(.05, .95)), data = samp)


y = df$val
# Plot the results
xx <- with(df, seq(min(y), max(y), len = 200))
phi.est <- logit(coef(fit2)[1], inverse = TRUE)
sd.est <- exp(coef(fit2)[3])
plot(xx, phi.est * dnorm(xx, Coef(fit2)[2], sd.est), col = "orange", type="l")
lines(xx, (1-phi.est) * dnorm(xx, Coef(fit2)[4], sd.est), col = "green")
abline(v = Coef(fit2)[c(2,4)], lty = 2, col = "orange")
## End(Not run)



fit2 <- vglm(val ~ 1, mix2normal(eq.sd = FALSE, qmu=c(.05, .95)), data = samp)


y = df$val
# Plot the results
xx <- with(df, seq(min(y), max(y), len = 200))
phi.est <- logit(coef(fit2)[1], inverse = TRUE)
sd.est1 <- exp(coef(fit2)[3])
sd.est2 <- exp(coef(fit2)[5])
plot(xx, phi.est * dnorm(xx, Coef(fit2)[2], sd.est1), col = "orange", type="l")
lines(xx, (1-phi.est) * dnorm(xx, Coef(fit2)[4], sd.est2), col = "green")
abline(v = Coef(fit2)[c(2)], lty = 2, col = "orange")
abline(v = Coef(fit2)[c(4)], lty = 2, col = "green")
## End(Not run)

d1 = density(df$val[df$roi ==0 ], n=75)
d2 = density(df$val[df$roi == 1], n=75)
lines(d1, col="orange", lwd=2)
lines(d2, col="green", lwd=2)


img.hist = hist(df$val, breaks=2000, 
 plot=FALSE)
 y = img.hist$counts
 x = img.hist$mids
 x = x[!is.na(y)];
 y = y[!is.na(y)]
 # 20 used for speed of example
 nawm_peak = get.last.mode(x, y, k=20, remove.tail=FALSE)

x = readNIfTI("100-362_20100126_1926_CT_2_CT_ROUTINE.nii.gz")
roi = readNIfTI("100-362_20100126_1926_CT_2_CT_ROUTINEROI.nii.gz")
df2 = data.frame(val = c(x))
df2$roi = c(roi)
sum(df2$roi[ df2$val < 20 | df2$val > 100])
df2 = df2[ df2$val >= 20 & df2$val <= 100, ]


backtoback = function(x, y, field= "counts"){
     ## using base
     h1 = hist(x, plot=FALSE)
     h2 = hist(y, plot=FALSE)
     all.field = c(h1[[field]], h2[[field]])
     h2[[field]] = - h2[[field]]
     hmax = max(all.field)
     # hmin = min(all.field)
     X = c(h1$breaks, h2$breaks)
     xmax = max(X)
     xmin = min(X)
     h1$counts = h1[[field]]
     h2$counts = h2[[field]]
     plot(h1, ylim=c(-hmax, hmax), col="green", xlim=c(xmin, xmax))
     lines(h2, col="blue")
}


img.hist = hist(df2$val, breaks=2000, 
 plot=FALSE)
 y = img.hist$counts
 x = img.hist$mids
 x = x[!is.na(y)];
 y = y[!is.na(y)]
 # 20 used for speed of example
 nawm_peak2 = get.last.mode(x, y, k=20, remove.tail=FALSE)


backtoback(df$val[df$roi == 1], df$val[df$roi == 0], field="density")
quartz()
backtoback(df2$val[df2$roi == 1], df2$val[df2$roi == 0], field="density")