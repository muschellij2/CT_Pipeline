#####################################
## Author: John Muschelli
## Date: January 20, 2014
## Purpose: Read in the AAL atlas labels and make R objects
## that can be used later for overlap metrics.
#####################################
#####################################
rm(list=ls())
library(R.matlab)
library(oro.nifti)
library(plyr)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
outdir = file.path(basedir, "results")
whichdir = "reoriented"


all.ids = list.dirs(basedir, recursive=FALSE, full.names=FALSE)
all.ids = all.ids[grepl("\\d\\d\\d-(\\d|)\\d\\d\\d", all.ids)]
redir = file.path(basedir, all.ids, whichdir)
redir = redir[file.exists(redir)]
nfiles = sapply(redir, function(x) 
  length(dir(path=x, pattern="ROI*.nii.*")))
redir = redir[nfiles > 0]
## those that have ROIs
all.ids = gsub(paste0(basedir, "/(.*)/", whichdir), "\\1", redir)
# all.ids = all.ids[file.exists(redir)]


spmdir = file.path(homedir, "spm8")

atlasdir = file.path(tempdir, "atlases")


tab.area = function(binimg, ind.list, keepall) {
  ## get overlap of indices
  raw.ind = which(binimg)
  raw.mat = sapply(ind.list, function(x) raw.ind %in% x)
  ## cs is sum of indices of overlap
  cs.raw = colSums(raw.mat)
  if (!keepall) cs.raw = cs.raw[cs.raw != 0]
  cs.raw = data.frame(nvox=cs.raw)  
  return(cs.raw)
}
bin.val = 0.95


get.pct = function(img, tempimg, Labels, ind.list,
  bin.val = 0.95,
  rescale=TRUE, 
  keepall=FALSE){
  ## you can pass either nifti iamge or filename 
  if (!inherits(img, 'nifti')){
    img = readNIfTI(img)
  }

  ### make sure same image
  at.pdim = pixdim(tempimg)[2:4]
  at.dim = dim(tempimg)

  dimg = dim(img)
  stopifnot(all( at.pdim == pixdim(img)[2:4]))
  stopifnot(all( at.dim == dim(img)))

  ## remove NAs and NaN - put to 0
  img[is.na(img)] = 0
  stopifnot(all(img >= 0))
  stopifnot(all(is.finite(img)))
  # img[is.nan(img)] = 0
  ###
  img.sum = sum(img)

  ## get binary image (any greater than 0)
  rawbin.img = img > 0
  n.vox = sum(rawbin.img)

  ### you may want to only call it a true voxel if greater than 
  ## bin.val (because of resampling)
  bin.img = img > bin.val
  nbin.vox = sum(bin.img)

  if (n.vox > 0){

    ## get overlap
    raw.tab = tab.area(rawbin.img, ind.list, keepall=keepall)
    bin.tab = tab.area(bin.img, ind.list, keepall= keepall)

  ## make sure everything went correctly
    stopifnot(sum(raw.tab$nvox) == n.vox)
    stopifnot(sum(bin.tab$nvox) == nbin.vox)

    img.ind = which(rawbin.img)
    ### these voxels have no categorization
    img.mat = sapply(ind.list, function(x) img.ind %in% x)
    # img.mat = img.mat[,rownames(raw.tab), drop=FALSE]
    ### weighted sum of area
    wsum = apply(img.mat, 2, function(inds){
      grab.ind = img.ind[inds]
      sum(img[grab.ind])
    })

    tol = 1e-5
    diff = abs(sum(wsum) - img.sum)
    stopifnot(diff < tol)
    wsum = data.frame(nvox=wsum)
    if (rescale){
      raw.tab = raw.tab   / n.vox
      bin.tab = bin.tab   / nbin.vox
      wsum  = wsum    / img.sum   
    }

    rownames(wsum) = rownames(raw.tab) = rownames(bin.tab) = Labels
  } else {
    wsum = data.frame(nvox=rep(0, length=length(Labels)+1))
    rownames(wsum) = c(Labels, "OutOfWindow")
    wsum$nvox[nrow(wsum)] = 1
    raw.tab = bin.tab = wsum
  }

  return(list(
    bin.val     = bin.val, 
    weighted.sum  = wsum, 
    raw.tab     = raw.tab,
    bin.tab     = bin.tab,
    n.vox       = n.vox,
    nbin.vox    =   nbin.vox,
    img.sum     =   img.sum))
}


collapse.res = function(res, add.binval=FALSE){
  cc = res[c("raw.tab", "bin.tab", "weighted.sum")]
  cc = lapply(cc, function(x){
    x$area = rownames(x)
    rownames(x) = NULL
    return(x)
  })
  c.res = merge(cc$raw.tab, cc$bin.tab, by="area", 
    all=TRUE, suffixes=c(".raw", ".bin"))
  c.res = merge(c.res, cc$weighted.sum, by="area", all=TRUE)
  colnames(c.res) = c("area", "raw", "bin", "weighted")
  nas = sapply(c.res, is.na)
  c.res[nas] = 0
  if (add.binval) c.res$bin.val = res$bin.val
  return(c.res) 
}




load(file.path(atlasdir, "All_FSL_Atlas_Labels.Rda"))
### contents
# tal.df, tal.img, 
#   mni.df, mni.img, 
#   hoxcort.df, hoxcort.img,
#   hoxsubcort.df, hoxsubcort.img, 




idir = 1
# for (idir in 1:2){
whichdir = "reoriented"

iid = all.ids[57]

mni.allres = tal.allres = NULL 

outfile = file.path(atlasdir, 
  paste0(whichdir, "_All_Atlas_ROI_Overlap_Measures.Rda"))

if (file.exists(outfile)) file.remove(outfile)

for (iid in all.ids){
  
### create directories
  redir = file.path(basedir, iid, whichdir)
  imgdir <- file.path(redir, "results")
  

  
  
  imgs = list.files(path=redir, 
                    full.names=TRUE, 
                    recursive=FALSE, 
                    pattern="^bws.*ROI.*nii")
  
  imgs = imgs[!grepl("affine9", imgs)]
  
  iimg = 1;
  x = imgs[iimg]
  # for (iimg in seq_along(imgs)){
  res.list = llply(imgs, function(x){
    # print(x)
    img.fname = x
    img = readNIfTI(img.fname)
    res = get.pct(img, tempimg=mni.img, Labels= mni.df$Label, 
      ind.list = mni.list, keepall=TRUE)
    mni.cres = collapse.res(res, add.binval=FALSE)
    mni.cres$fname = x

    Labels = apply(tal.df[, 
      c("Area", "Lobe", "Lobule", "Tissue_Type", "Broadmann")], 
      1, paste, sep='', collapse=".")
    res = get.pct(img, tempimg=tal.img, Labels= Labels, 
      ind.list = tal.list, keepall=TRUE)
    tal.cres = collapse.res(res, add.binval=FALSE)
    tal.cres$fname = x

    # iimg <<- iimg + 1
    return(list(mni.cres=mni.cres, tal.cres = tal.cres))
  }, .progress= "text")
  
  mni.res = lapply(res.list, function(x) x$mni.cres)
  tal.res = lapply(res.list, function(x) x$tal.cres)


  mni.res = do.call("rbind", mni.res)
  mni.res$fname = basename(mni.res$fname)
  mni.res$fname = gsub("\\.gz$", "", mni.res$fname)
  mni.res$fname = gsub("\\.nii$", "", mni.res$fname)
  mni.res$fname = gsub("^2mm_", "", mni.res$fname)
  mni.res$fname = gsub("affine(9|12)_", "", mni.res$fname)

  tal.res = do.call("rbind", tal.res)
  tal.res$fname = basename(tal.res$fname)
  tal.res$fname = gsub("\\.gz$", "", tal.res$fname)
  tal.res$fname = gsub("\\.nii$", "", tal.res$fname)
  tal.res$fname = gsub("^2mm_", "", tal.res$fname)
  tal.res$fname = gsub("affine(9|12)_", "", tal.res$fname)

  mni.res = mni.res[ mni.res$raw > 0, ]
  tal.res = tal.res[ tal.res$raw > 0, ]

  tal.allres = rbind(tal.allres, tal.res)
  mni.allres = rbind(mni.allres, mni.res)
  print(iid)
}

  save(mni.allres, tal.allres, 
       file=outfile)
# }
# cres = cres[order(cres$weighted),]
# }