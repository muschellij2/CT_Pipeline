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

jhut2.allres = jhut1.allres = mni.allres = tal.allres = NULL 

outfile = file.path(atlasdir, 
  paste0(whichdir, "_All_Atlas_ROI_Overlap_Measures.Rda"))

if (file.exists(outfile)) file.remove(outfile)

for (iiid in seq_along(all.ids)){
  

  iid = all.ids[iiid]
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
    my.pb = txtProgressBar(min = 0, max=4, style=3)
    
    img.fname = x
    img = readNIfTI(img.fname)
    res = get.pct(img, tempimg=mni.img, Labels= mni.df$Label, 
      ind.list = mni.list, keepall=TRUE)
    mni.cres = collapse.res(res, add.binval=FALSE)
    mni.cres$fname = x

    ### talairach
    setTxtProgressBar(my.pb, 1)
    Labels = apply(tal.df[, 
      c("Area", "Lobe", "Lobule", "Tissue_Type", "Broadmann")], 
      1, paste, sep='', collapse=".")
    res = get.pct(img, tempimg=tal.img, Labels= Labels, 
      ind.list = tal.list, keepall=TRUE)
    tal.cres = collapse.res(res, add.binval=FALSE)
    tal.cres$fname = x

    setTxtProgressBar(my.pb, 2)
    ## EVE Template 1
    res = get.pct(img, tempimg=jhut1.img, Labels= jhut1.df$Label, 
      ind.list = jhut1.list, keepall=TRUE)
    jhut1.cres = collapse.res(res, add.binval=FALSE)
    jhut1.cres$fname = x    

    setTxtProgressBar(my.pb, 3)
    ## EVE Template 2
    res = get.pct(img, tempimg=jhut2.img, Labels= jhut2.df$Label, 
      ind.list = jhut2.list, keepall=TRUE)
    jhut2.cres = collapse.res(res, add.binval=FALSE)
    jhut2.cres$fname = x        
    setTxtProgressBar(my.pb, 4)
    close(my.pb)
    # iimg <<- iimg + 1
    return(list(mni.cres=mni.cres, 
      tal.cres = tal.cres,
      jhut1.cres = jhut1.cres,
      jhut2.cres = jhut2.cres
      ))
  }, .progress= "text")
  
  mni.res = lapply(res.list, function(x) x$mni.cres)
  tal.res = lapply(res.list, function(x) x$tal.cres)
  jhut1.res = lapply(res.list, function(x) x$jhut1.cres)
  jhut2.res = lapply(res.list, function(x) x$jhut2.cres)

  make.mat = function(res){
    res = do.call("rbind", res)
    res$fname = basename(res$fname)
    res$fname = gsub("\\.gz$", "", res$fname)
    res$fname = gsub("\\.nii$", "", res$fname)
    res$fname = gsub("^2mm_", "", res$fname)
    res$fname = gsub("affine(9|12)_", "", res$fname)
    res = res[ res$raw > 0, ]
    return(res)
  }

  mni.res = make.mat(mni.res)
  tal.res = make.mat(tal.res)
  jhut1.res = make.mat(jhut1.res)
  jhut2.res = make.mat(jhut2.res)


  tal.allres = rbind(tal.allres, tal.res)
  mni.allres = rbind(mni.allres, mni.res)
  jhut1.allres = rbind(jhut1.allres, jhut1.res)
  jhut2.allres = rbind(jhut2.allres, jhut2.res)  
  print(iid)
}

  save(mni.allres, tal.allres, 
      jhut1.allres, jhut2.allres, 
       file=outfile)

# }
# cres = cres[order(cres$weighted),]
# }