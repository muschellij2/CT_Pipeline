############################################
## Get density of all HU
############################################
############################################
rm(list=ls())
library(plyr)
library(cttools)
library(fslr)
library(mgcv)
library(methods)
homedir = "/Applications"
rootdir = file.path("/Volumes/DATA_LOCAL",
    "Image_Processing")
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = file.path("/legacy/dexter/disk2/smart", 
    "stroke_ct", "ident")
}
progdir = file.path(rootdir, "programs")
segdir = file.path(progdir, "Segmentation")

basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")
outdir = file.path(basedir, "results")

#### load voxel data

outfile = file.path(outdir, 
    "Reseg_111_Filenames.Rda")
xxx = load(file = outfile)

iimg <- suppressWarnings({
    as.numeric(Sys.getenv("SGE_TASK_ID"))
    })
if (is.na(iimg)) iimg = 34
## 2,12,17,34,37, 46, 48, 68, 
## 70,81, 85, 86, 87, 99
#  has variable slice thickness
## 15 is not correct for sthickness
## 17 & 87 worst - has overlapping 
## slice somewhat
## 75 has all negative
## 71 has no position data
## 13,71,101 has spacing



# fdf = fdf[c(2,12,17,34,37, 
    # 46, 48, 68, 70,81, 
    #85, 86, 87, 99),]
# dcmtables[, '0018-1152-Exposure']
N = nrow(fdf)
ds = vector(mode = "list", length = N)
vs = ds
for (iimg in seq(N)){
    
    runx = x = fdf[iimg,]
    sortdir = file.path(x$iddir, "Sorted")

    # run_model = function(x, fpr.stop = .1){
    fname = xfname = nii.stub(x$img, bn=TRUE)
    print(iimg)
    fname = file.path(x$outdir, 
        paste0("Reseg_", xfname, 
            "_hu.Rda"))
    load(fname)
    vs[[iimg]] = vals
    ds[[iimg]] = dvals
    # print(warnings())
}

names(vs) = names(ds) = fdf$id

ltrain = fdf$group %in% "Train"
ltest = !ltrain


############################
# Get ranges for plotting
############################
max_y = sapply(ds, function(x) {
    range(x$y)
})
ry = range(max_y)

# plot the densities
# plot(ds[[1]], xlim = c(0, 100), ylim = ry)
# for (i in 2:N) {
#     lines(ds[[i]])
# }
# abline(v = c(40, 80), col = "red")

############################
# Collapse all values
all_vals = unlist(vs)
d = density(all_vals)
cdf = ecdf(all_vals)
plot(d)
abline(v = c(40, 80), col = "red")


train_vals = unlist(vs[ltrain])
dtrain = density(train_vals)
train_cdf = ecdf(train_vals)
plot(dtrain)
abline(v = c(40, 80), col = "red")


scaled = lapply(ds, function(x) {
    x$y = x$y / max(x$y)
    return(x)
})
plot(scaled[[1]])
for (i in 2:N) {
    lines(scaled[[i]])
}