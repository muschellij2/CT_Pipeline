rm(list=ls())
library(cttools)
library(fslr)
library(plyr)
homedir = "~/"
rootdir = "~/CT_Registration/"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
homedir = file.path(rootdir, "Rorden_data")
progdir = file.path(rootdir, "programs")
basedir = file.path(rootdir, "Registration")

# scandir = file.path(homedir, "ct_scans")
regdir = file.path(homedir, "registered_ct_scans")
# setwd(homedir)
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")
outdir = file.path(basedir, "results")

template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)
dtemp = dim(temp)


maskfile = file.path(tempdir, "scct_mask.nii.gz")
tempmask = readNIfTI(maskfile)

files = list.files(regdir, pattern = ".*SS_.*.nii.gz", 
	full.names=TRUE)
ifile = 1

pb = txtProgressBar(min=0, max=length(files), style=3)
i = 0
imgs = lapply(files, function(x) {
    i <<- i + 1
    setTxtProgressBar(pb, i)
	readNIfTI(x, reorient=FALSE)
})
close(pb)

mat = sapply(imgs, c)
out = which(mat > 100, arr.ind=TRUE)
mat[out] = 0

mat[mat == 0] = NA

rmean = rowMeans(mat, na.rm=TRUE)
rsd = rowSds(mat, na.rm=TRUE)

Ns = rowSums(!is.na(mat))

sd.img = newnii(temp)
sd.img@.Data= array(rsd, dim = dtemp)
sd.img = cal_img(sd.img)


mn.img = newnii(temp)
mn.img@.Data= array(rmean, dim = dtemp)
mn.img = cal_img(mn.img)

N.img = newnii(temp)
N.img@.Data= array(Ns, dim = dtemp)
N.img = cal_img(N.img)


# n = length(imgs)
ch_rsd = sqrt(rsd^2 * (Ns+1) / Ns)
ch_sd.img = newnii(temp)
ch_sd.img@.Data= array(ch_rsd, dim = dtemp)
ch_sd.img = cal_img(ch_sd.img)



fname = file.path(regdir, "Mean_Image")
writeNIfTI(mn.img, fname)

fname = file.path(regdir, "SD_Image")
writeNIfTI(sd.img, fname)

fname = file.path(regdir, "CH_SD_Image")
writeNIfTI(ch_sd.img, fname)

fname = file.path(regdir, "N_Image")
writeNIfTI(N.img, fname)

