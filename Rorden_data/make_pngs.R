rm(list=ls())
library(cttools)
library(fslr)
library(plyr)
library(extrantsr)
homedir = "~/"
rootdir = "~/CT_Registration/"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
rootdir = path.expand(rootdir)
homedir = file.path(rootdir, "Rorden_data")
progdir = file.path(rootdir, "programs")
basedir = file.path(rootdir, "Registration")

# scandir = file.path(homedir, "ct_scans")
regdir = file.path(homedir, "registered_ct_scans")
scandir = file.path(homedir, "ct_scans")
outdir = file.path(homedir, "pngs")
avgdir = file.path(homedir, "Averages")
# setwd(homedir)
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")
# outdir = file.path(homedir, "SyN_Registered_Scans")

template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)
dtemp = dim(temp)

ss_template = file.path(tempdir, "scct_unsmooth_SS_0.01.nii.gz")
sstemp = readNIfTI(ss_template)

pngname = paste0("Template_Slice.png")
pngname = file.path(outdir, pngname)
ss = dropEmptyImageDimensions(sstemp)
png(pngname, res = 600, units ="in", height=7, width = 7, type ="cairo")
  image(ss, z = ceiling(nsli(ss)/2), plot.type = "single")
dev.off()

mean.fname = file.path(avgdir, "Rorden_Mean_Image.nii.gz")
meanimg = readNIfTI(mean.fname)
pngname = paste0("Mean_Slice.png")
pngname = file.path(outdir, pngname)
ss = dropEmptyImageDimensions(meanimg)
png(pngname, res = 600, units ="in", height=7, width = 7, type ="cairo")
image(ss, z = ceiling(nsli(ss)/2), plot.type = "single")
dev.off()

sd.fname = file.path(avgdir, "Rorden_SD_Image.nii.gz")
sdimg = readNIfTI(sd.fname)
pngname = paste0("SD_Slice.png")
pngname = file.path(outdir, pngname)
ss = dropEmptyImageDimensions(sdimg)
png(pngname, res = 600, units ="in", height=7, width = 7, type ="cairo")
image(ss, z = ceiling(nsli(ss)/2), plot.type = "single")
dev.off()

ss_files = list.files(path=scandir, full.names=TRUE, 
                      pattern="^.*_SS_0.01\\.nii")

ss_files =  ss_files[1:3]
iname = 1

for (iname in seq_along(ss_files)){
  pngname = paste0("Reg_Image", iname, ".png")
  pngname = file.path(outdir, pngname)
  png(pngname, res = 600, units ="in", height=7, width = 7, type ="cairo")
    ss = readnii(ss_files[iname])
    ss = dropEmptyImageDimensions(ss)
    image(ss, z = ceiling(nsli(ss)/2), plot.type = "single")
  dev.off()
}

ifile = 1
