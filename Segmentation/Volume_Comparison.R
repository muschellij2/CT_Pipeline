##################################################################
## This code is for plotting volumes
##
## Author: John Muschelli
## Last updated: May 20, 2014
##################################################################
##################################################################
rm(list=ls())
library(ggplot2)

lm_eqn = function(df){
    m = lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(coef(m)[1], digits = 3), 
              b = format(coef(m)[2], digits = 3), 
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
}

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

outfile = file.path(outdir, 
    paste0("Cutoff_predictions.Rda"))
load(file = outfile)

types = c("_zval2", '_zval2_medztemp')
cutter = c("", "_dice")

vals = c(outer(types, cutter, paste0))
vals = c(outer(vals, c("_pred"), paste0))
vals = paste0("vol", vals)
xvals = vals

vals = paste0("vol", types, "_raw")
vals = c(vals, xvals)

filename = file.path(outdir, 
    paste0("Result_Formats.Rda"))
load(filename)

nopred = run.ind
ffdf = fdf[-nopred, ]
nr = nrow(ffdf)
valid.ind = ceiling(nr/2)
test.ind = seq( valid.ind +1, nr)
valid.ind = seq(1, valid.ind)

group = "Training"
if (group == "Training"){
    subset.ind = valid.ind
}
if (group == "Test"){
    subset.ind = test.ind
}

vols = ffdf[, c("id", vals, "truevol")] 

dat = vols[subset.ind,]
ival = vals[1]
dat$y = dat$truevol
pdfname = file.path(outdir, "Cutoff_Volume_Comparison.pdf")
pdf(pdfname)
for (ival in vals){
    dat$x = dat[, ival]
    lab = gsub("vol_", "", ival)
    ztemp = grepl("ztemp", ival)
    raw = grepl("raw", ival)
    p = ggplot(dat, aes(y=y, x=x)) 

    p = p + geom_point() + 
        geom_abline(aes(intercept=0, slope = 1, 
        colour="X=Y"))
    p = p + geom_smooth(aes(colour="Linear Fit"), 
        method="lm",  se=FALSE)

    p = p + ylab("Manual Volume (mL)") + 
        xlab("Predicted Volume (mL)")
    # p = p + ggtitle(ival)
    lab = paste0("Adding", ifelse(raw, " Raw", " Thresholded"), 
        " Probabilities with General", 
        ifelse(ztemp, '+Subject-Specific', ""), " Voxel Selection")
    p = p + ggtitle(lab)
    # p = p + geom_smooth(col="blue", se=FALSE)
    p1 = p + geom_text(aes(x = 35, y = 100, 
        label = lm_eqn(dat)), parse = TRUE)

    tsize = 16
    p1 = p1 + theme(legend.position = c(0.25, 0.65),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill="transparent", 
                                  color="transparent"),
        legend.text = element_text(size=tsize+2), 
        legend.title = element_text(size=tsize),
        title = element_text(size=tsize),
        strip.text = element_text(size = tsize+4),
        axis.text  = element_text(size=tsize-2))  +
  scale_colour_manual("", 
                      values = c("X=Y"="black",
                        "Linear Fit" ="red"),
                      guide = guide_legend(reverse=TRUE))   
    # plot(truevol ~ ., data=dat[c("truevol", ivol)] )
    # mod = lm(truevol ~ ., data=dat[c("truevol", ivol)])
    # abline(a=0, b=1)
    # abline(mod, col="red")
    print(p1)
    print(ival)
}
dev.off()

dat$y = dat$x = NULL
    # print(iscen)
# }

