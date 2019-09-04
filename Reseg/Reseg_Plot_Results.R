####################################
## This code is for aggregate models
## CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
####################################
#####################################
rm(list=ls())
library(methods)
library(plyr)
library(dplyr)
library(cttools)
library(fslr)
library(ggplot2)
library(reshape2)
set.seed(20150518)
homedir = "/Applications"
rootdir = file.path("/Volumes/DATA_LOCAL", 
    "Image_Processing")
if (Sys.info()[["user"]] %in% 
    "jmuschel") {
  homedir = "~"
  rootdir = file.path(
    "/legacy/dexter/disk2/smart", 
    "stroke_ct", "ident")
}
progdir = file.path(rootdir, 
    "programs")
basedir = file.path(rootdir, 
    "Registration")
tempdir = file.path(rootdir, 
    "Template")
atlasdir = file.path(tempdir, 
    "atlases")
outdir = file.path(basedir, 
    "results")

segdir = file.path(progdir, 
    "Reseg")
source(file.path(segdir, 
    "Reseg_performance_functions.R"))
source(file.path(segdir, 
    "ggplot_smooth_func.R"))

mypng = function(pngname, 
    height = 7, width =7){
    png(pngname,
        height = height,
        width = width,
        units = "in",
        type = "cairo", 
        res = 600)
}


tsize = 16
myopts = theme(legend.position = c(.36, .65),
        legend.key = element_rect(fill="transparent", 
                                  color="transparent"),
        legend.text = element_text(size=tsize+4), 
        legend.title = element_text(size=tsize),
        title = element_text(size=tsize),
        plot.title = element_text(hjust = 0.8),        
        strip.text = element_text(size = tsize+4),
        axis.text  = element_text(size=tsize-2))

# 



outrda = file.path(outdir, 
    paste0("Reseg_Results.Rda"))
xxx = load(outrda)

xlong = long

long = filter(long, 
    cutoff %in% c("cc", "scc"))
long$cutoff = revalue(long$cutoff, 
    c("cc"= "Unsmoothed",
    "scc" = "Smoothed")
    )
long = mutate(long, 
    mean = (tvol + evol) /2,
    diff = tvol - evol
    )
slong = filter(long, 
    cutoff %in% c("Smoothed"),
    group %in% c("Test", "Validation"))
relev = c("logistic" = "Logistic",
        "lasso" = "LASSO",
        "gam" = "GAM",
        "rf" = "Random\nForest")
slong$mod = revalue(slong$mod, 
    relev
    )
slong$mod = factor(slong$mod,
    levels = relev)

nlong = filter(slong, app %in% "Native")
llong = select(nlong, mod, 
    dice, sens, accur,
    spec, iimg, group)
llong = melt(llong, 
    id.vars = c("iimg", "group", "mod"))
relev2 = c("dice" = "Dice Similarity Index",
        "accur" = "Accuracy",
        "sens" = "Sensitivity",
        "spec" = "Specificity")
llong$variable = revalue(llong$variable, 
    relev2
    )
llong$variable = factor(llong$variable, 
    levels = relev2)

baplot = nlong %>% 
    ggplot(aes(x = mean, y = diff)) + 
        geom_point() + facet_wrap(~ mod) +
        geom_smooth(se = FALSE) +
        geom_hline(yintercept = 0,
            col = "red")
baplot = baplot + myopts
baplot = baplot + 
    xlab("ICH Volume Mean (Manual + Estimated)/2") +
    ylab("ICH Volume Difference (Manual - Estimated)") 

pngname = file.path(outdir, 
    "Reseg_BA_Plots.png")
mypng(pngname)
    print(baplot)
dev.off()    



# f = long %>% filter(
#         mod %in% c("rf", "logistic")) %>%
#     select(dice, mod, app, 
#         cutoff, iimg, group)
# dwide = dcast(f, 
#     iimg + app + group + cutoff  ~ mod,
#     value.var = "dice")
# dwide = mutate(dwide, 
#     mean = (logistic + rf)/2,
#     diff = (logistic - rf))
# dtest = filter(dwide, group %in% "Validation")

# scat = ggplot(aes(x = logistic, y = rf), 
#     data = dtest) + 
#     geom_point() + 
#     facet_wrap( ~ app + cutoff)
# scat

# ba = ggplot(aes(x = mean, y = diff), 
#     data = dwide) + 
#     geom_point() + 
#     facet_grid(app ~ cutoff) + 
#     geom_hline(aes(y = 0)) + 
#     geom_smooth(se = FALSE)
# ba 
# ba %+% filter(dwide, mean > 0.2)


# g = ggplot(aes(x = mod, colour = cutoff), 
#     data = long[ 
#         long$group == "Test",]) + 
#     geom_boxplot() + facet_wrap( ~ app)
# gdice = g + aes(y = dice)

# print(gdice)


diffs = nlong  %>% 
    group_by(mod) %>% 
    summarise(
        mn_diff = mean(diff),
        med_diff = median(diff),
        sd_diff = sd(diff),
        se_diff = sd(diff)/n(),
        mn_rat = mean(evol/tvol),
        med_rat = median(evol/tvol),
        sd_rat = sd(evol/tvol),
        se_rat = sd(evol/tvol)/n(),
        mn_pct_diff = mean(diff/tvol),
        med_pct_diff = median(diff/tvol),
        sd_pct_diff = sd(diff/tvol),
        se_pct_diff = sd(diff/tvol)/n()
        ) %>% 
    as.data.frame

outdiff = file.path(outdir, 
    paste0("Reseg_Differences_Results.Rda"))
save(diffs, file = outdiff)

g = ggplot(aes(x = mod), 
    data = nlong) + 
    geom_boxplot() +
    myopts 
ogdice = g + aes(y = dice)
ogdice = ogdice + xlab("Model") + 
    ylab("Dice Similarity Index")
pngname = file.path(outdir, 
    "Reseg_Dice_Comparison.png")
mypng(pngname)
    print(ogdice)
print(ogdice)
dev.off()    

# g = ggplot(aes(x = mod, colour = app), 
#     data = slong) + 
#     geom_boxplot() 
# ogdice2 = g + aes(y = dice)
# print(ogdice2)

# g = ggplot(aes(x = mod, colour = group), 
#     data = slong) + 
#     geom_boxplot() + facet_wrap( ~ app)
# gdice = g + aes(y = dice)
# print(gdice)

g = ggplot(aes(x = mod, colour = group), 
    data = nlong) + 
    geom_boxplot() 
gdice = g + aes(y = dice)
# print(gdice)


g = ggplot(aes(x = variable,
    y = value, colour = mod), 
    data = llong) + 
    geom_boxplot()
# print(g)

g = ggplot(aes(x = mod,
    y = value), 
    data = llong) + 
    geom_boxplot() 
gd = g + facet_wrap(~ variable)
gd2 = g + facet_wrap(~ variable, 
    scales = "free_y")
# print(gd)
# print(gd2)

dsens = filter(llong,
    variable %in% c("Sensitivity", 
        "Dice Similarity Index"))
g = ggplot(aes(x = mod,
    y = value), 
    data = dsens) + 
    geom_boxplot() 
gd = g + facet_wrap(~ variable)
# print(gd)

# mean_vol = floor(mean(long$evol))
# long$mevol = long$evol - mean_vol
getfunc = function(x){
    rmse = sqrt(mean((x$tvol - x$evol)^2))
    pct_diff = (x$tvol - x$evol)/x$tvol
    pct_rmse = sqrt(mean(pct_diff^2))

    mod = lm(tvol ~ evol, data =x)
    smod = summary(mod)
    cc = coef(mod)
    volcor = unlist(
        summarise(x, cor(tvol, evol))
        )
    volcor2 = unlist(
        summarise(x, cor(tvol, evol,
            method = "spearman"))
        )        
    ct = cor.test(x$tvol, x$evol)
    names(volcor) = "cor"
    names(volcor2) = "scor"
    names(cc) = c("Intercept", "Slope")
    se = c(coef(smod)[, "Std. Error"])
    names(se) = c("se.Intercept", 
        "se.Slope")
    ci = confint(mod)
    colnames(ci) = c("lower", "upper")
    ci = unlist(data.frame(ci))
    names(ci) = gsub("1$", ".Intercept", 
        names(ci))
    names(ci) = gsub("2$", ".Slope", 
        names(ci))  
    ct.ci = ct$conf.int      
    names(ct.ci) = c("cor.lower", "cor.upper")
    cc = c(cc, 
        slope_diff = abs(cc["Slope"]-1),
        se,
        ci,
        volcor, volcor2,
        rmse = rmse,
        pct_rmse = pct_rmse,
        ct.ci
        )
    return(c(cc))
}
native = filter(slong, app %in% "Native")
gcorrs = ddply(native, .(group, mod), 
    getfunc)
corrs = ddply(native, .(mod), getfunc)
rownames(corrs) = names(relev)

outcorr = file.path(outdir, 
    paste0("Reseg_Correlation_Results.Rda"))
save(corrs, file = outcorr)




rd = range(c(native$evol, native$tvol))
g = ggplot(aes(y = tvol, x = evol), 
    data = native) + 
    geom_point() + theme_bw() +myopts
g = g +
 geom_segment(
    aes(x = 0, y = 0,
        xend = rd[2], yend = rd[2],
    colour = "X = Y Line"),
  size = .75)
greg = g
g = g + stat_smooth_func(geom="text",
    method="lm",
    hjust=0,parse=TRUE,
    colour = "blue",
    ymultiplier = c(0.9, 0.8),
    n = 10000) +
  geom_smooth(aes(colour = "Linear Fit"),
    method="lm", se=FALSE,
    size = .75) 
g = g + 
scale_colour_manual("", 
values = c("X = Y Line" ="deeppink", 
    "Linear Fit"="blue"),
guide = guide_legend(reverse=TRUE)) +
theme(legend.position = c(.36, .74))

g = g + theme(
legend.background = 
element_rect(fill="transparent"))
g = g + xlab("Automatic Segmentation Volume (cc)") +
    ylab("Manual Segmentation Volume (cc)")
g = g + 
    scale_x_continuous(
        breaks = seq(0, 150, by = 30),
        limits = c(0, 150)
        ) +
    scale_y_continuous(
        breaks = seq(0, 150, by = 30),
        limits = c(0, 150)        
        )
g2 = g + facet_grid(mod ~ group)
gg = g
g = g + facet_wrap(~ mod)

d = data.frame(
    mod = unique(native$mod)
    )
levels(d$mod) = levels(native$mod)
d = d[ order(d$mod),,drop = FALSE]
d$label = LETTERS[seq(nrow(d))]

gnolet = g
g = g + geom_text(
    data=d, x = 120, y = 15, 
    size=20, aes(label=label), 
    colour="black")

pngname = file.path(outdir, 
    "Reseg_Volume_Comparison.png")
mypng(pngname)
    print(g)
dev.off()


pngname = file.path(outdir, 
    "Reseg_Volume_Comparison_noletters.png")
mypng(pngname)
    print(gnolet)
dev.off()



pngname = file.path(outdir, 
    "Reseg_Volume_Logistic.png")
mypng(pngname)
print({
    g %+% filter(native, mod %in% "Logistic")
    })
dev.off()



pngname = file.path(outdir, 
    "Reseg_Volume_Comparison_Long.png")
mypng(pngname,
    height = 3.5, 
    width = 3.5*4)
gg = g + facet_wrap(~mod, nrow = 1)
gg + theme(legend.position = c(0.42, 0.17))
dev.off()


print(g2)

glong = filter(long, group == "Test")
rigid = filter(glong, app %in% "Rigid")
native = filter(glong, app %in% "Native")


vlong = filter(long, group == "Validation")
rigid_valid = filter(vlong,
    app %in% "Rigid")
native_valid = filter(vlong,
    app %in% "Native")



g = ggplot(aes(y = tvol, x = evol), 
    data = rigid) + 
    geom_point() + 
    facet_wrap(~ cutoff + mod)
g = g + stat_smooth_func(geom="text",
    method="lm",
    hjust=0,parse=TRUE) +
  geom_smooth(method="lm",se=FALSE) +
  geom_abline(aes(intercept = 0, 
        slope = 1))
g
g %+% native
g %+% native_valid
g %+% rigid_valid


ba = ggplot(aes(y = tvol, x = evol), 
    data = rigid) + 
    geom_point() + 
    facet_wrap(~ cutoff + mod)
ba = ba + stat_smooth_func(geom="text",
    method="lm",
    hjust=0,parse=TRUE) +
  geom_smooth(method="lm",se=FALSE) +
  geom_abline(aes(intercept = 0, 
        slope = 1))
g
g %+% native
g %+% native_valid
g %+% rigid_valid


######################################
# Random forest dices
######################################
rf = filter(long, mod == "rf")
rf = select(rf, iimg, cutoff, dice, group, 
    app)
rf = mutate(rf, 
    smooth = c("Unsmoothed", "Smoothed")[
    grepl("^s", cutoff) + 1]
    )
rf$cutoff = gsub("^s", "", rf$cutoff)

rf1 = reshape(rf, 
    direction = "wide",
    idvar = c("iimg", "group", 
        "smooth", "app"),
    timevar = "cutoff"
    )
colnames(rf1) = gsub("dice[.]", "",
    colnames(rf1))

rf2 = reshape(rf, 
    direction = "wide",
    idvar = c("iimg", "group", 
        "cutoff", "app"),
    timevar = "smooth"
    )
colnames(rf2) = gsub("dice[.]", "",
    colnames(rf2))

# gg2 = ggplot(
#     aes(x = Smoothed, y = Unsmoothed, 
#         colour = group), 
#     data = rf2) + 
#     geom_point() + 
#     facet_wrap(~ cutoff + app) + 
#     geom_abline(aes(intercept = 0, 
#         slope = 1))
# gg2

# gg = ggplot(
#     aes(x = pred, y = cc, colour = group), 
#     data = rf1) + 
#     geom_point() + 
#     facet_wrap(~ smooth + app) + 
#     geom_abline(aes(intercept = 0, 
#         slope = 1))
# gg
