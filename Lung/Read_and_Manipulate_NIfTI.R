rm(list=ls())
library(cttools)
library(fslr)
library(plyr)
options(matlab.path='/Applications/MATLAB_R2013b.app/bin')

rootdir = "~/Dropbox/_LoriLungs/NIfTI"
rootdir = path.expand(rootdir)

basedir = rootdir

niis = list.files(path=basedir, pattern="*.nii.gz", full.names=TRUE)

nii = niis[2]

threshdir = file.path(basedir, "thresholds")
if (!file.exists(threshdir)){
    dir.create(threshdir)
}

# img = readNIfTI(nii, reorient=FALSE)
img = readNIfTI(nii)

########
newnifti = function(img, arr){
    x = img
    x@.Data = arr
    x = cal_img(x)
}

########## cal_img is from fslr
######## will reset cal_max and cal_min slots
# mask = newnifti(img, img > 100 & img < 1000)
# orthographic(img)
# orthographic(mask)
# mask.overlay(img, mask)

thresh = seq(-900, -700, by=25)
stub = basename(nii.stub(nii))

pdfname = file.path(threshdir, paste0(stub, "_thresholds.pdf"))
pdf(pdfname)
for (ithresh in thresh){
    orthographic(cal_img(img > ithresh), 
        text=paste0("Threshold is ", ithresh),
        useRaster = TRUE)
}
dev.off()

mask = newnifti(img, img > -800 & img <= 100)
orthographic(mask)

orthographic(cal_img(img > 100))