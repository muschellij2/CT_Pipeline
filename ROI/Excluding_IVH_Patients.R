#################################
# Regressions with % of ROI
# Author: John Muschelli
#################################
rm(list=ls())
library(cttools)
library(scales)
library(RColorBrewer)
library(data.table)
library(ggplot2)
library(grid)
library(plyr)
library(fslr)
homedir = "/Applications"
rootdir = "~/CT_Registration"
basedir = file.path(rootdir, "data")
outdir = basedir
figdir = file.path(rootdir, "CT_Pipeline", "figure")
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
  basedir = file.path(rootdir, "Registration")
  outdir = file.path(basedir, "results")
}
progdir = file.path(rootdir, "programs")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")


img_cut = function(img, breaks, ...){
  cuts = cut(img, breaks=breaks, ...)
  # cuts = factor(cuts, levels)
  levs = levels(cuts)
  cuts = as.numeric(cuts)
  # res.p[ rs > ncut ] = cuts
  img = niftiarr(img, array(cuts, dim=dim(img)))
  return(list(img=img, levs=levs))
}

lr_symm = function(img){
  dimg = dim(img)
  max.slice = dimg[1]  
  mid.slice = (max.slice+1)/2
  
  w = which(img > 0, arr.ind=TRUE)
  ## 20 - then 160, 90 - 20 + 90
  ## if 160 then 90 - 160 + 90
  w[, 1] = 2 * mid.slice - w[,1]
  w = w[ w[, 1] > 0 & w[, 1] < max.slice, ]
  img[w] = 1
  
  img = (img > 0)*1
  img = newnii(img)
}

# allres = allres
make.pvalimg = function(pvalimg, runlist = NULL){
  pvalimg.tab = llply(runlist, function(x) {
    x = area_pct(pvalimg, ind.list=x, keepall=TRUE)		
    x$nvox = x$nvox/sum(x$nvox) * 100
    x$roi_pct_any = x$roi_pct_any * 100
    x$roi_mean_pct = x$roi_mean_pct * 100		
    x = x[order(x$nvox, decreasing=TRUE), , drop=FALSE]
    x$area = rownames(x)
    x
  }, .progress= "text")
  
  names(pvalimg.tab) = names(runlist)
  return(pvalimg.tab)
}



atfile = file.path(atlasdir, "All_FSL_Atlas_Labels.Rda")

x = load(file=atfile)


lists = list(mni.list, jhut1.list, jhut2.list)
names(lists) = c("MNI", "EVE_1", "EVE_2")
sublists = list(jhut1.list, jhut2.list)


sublists = lapply(sublists, function(x) {
  area = names(x)
  x[grep("GLOBUS_PALLIDUS|THALAMUS|PUTAMEN", area)]
})

sublists = lapply(sublists, function(x) {
  xx = unlist(x)
  area = names(xx)
  area = gsub("_left\\d*", "", area)
  area = gsub("_right\\d*", "", area)
  uarea = unique(area)
  x = lapply(uarea, function(aname){
    ind = which(area %in% aname)
    xx[ind]
  })
  names(x) = uarea
  x
})


col.lists = list(jhut1.list, jhut2.list)
names(col.lists) = c("EVE_1", "EVE_2")

col.lists = lapply(col.lists, function(x) {
  area = names(x)
  area = gsub("_left", "", area)
  area = gsub("_right", "", area)
  uarea = unique(area)
  res = lapply(uarea, function(aname){
    ind = which(area %in% aname)
    xx = sort(unlist(x[ind]))
    names(xx) = NULL
    # print(ind)
    xx
    # xx[ind]
  })
  names(res) = uarea
  res
})

rm(list=x)


area_pct = function(img, ind.list, keepall) {
  ## get overlap of indices
  raw.mat = sapply(ind.list, function(x) sum(img[x]))
  any.mat = sapply(ind.list, function(x) mean(img[x] > 0))
  mn.mat = sapply(ind.list, function(x) mean(img[x]))
  names(raw.mat) = names(ind.list)
  ## cs is sum of indices of overlap
  cs.raw = data.frame(nvox=raw.mat, roi_pct_any = any.mat,
                      roi_mean_pct = mn.mat) 
  rownames(cs.raw) = names(ind.list)
  if (!keepall) cs.raw = cs.raw[rowSums(cs.raw) != 0, , 
    drop=FALSE]
  return(cs.raw)
}

eve = col.lists[[2]]
ind = eve[[which(names(eve) == "Background")]]

whichdir = "reoriented"
outcome = "NIHSS"


get.id = function(x){
  ss = strsplit(x, "_")
  ss = sapply(ss, head, 1)
  ss = gsub(".*(\\d\\d\\d-.*)", "\\1", ss)
  ss
}

id_to_pname = function(x){
  as.numeric(gsub("-", "", x))
}

nkeeps = c(1000, 2000, 3000, .001, 0.01, 0.05)

demog = read.csv(file=file.path(basedir, "Demog_NIHSS_Mask.csv"), 
                 stringsAsFactors=FALSE)
demog$Base_ICH_10 = demog$Diagnostic_ICH /10

demog$Clot_Location_RC = gsub("Palidus", "Pallidus", 
                              demog$Clot_Location_RC )
demog$Clot_Location_RC = factor(demog$Clot_Location_RC, 
                                levels= 
                                  c("Lobar", "Globus Pallidus", 
                                    "Putamen", "Thalamus"))
demog$LOC = demog$Clot_Location_RC


all.demog = demog

###############################################
# Load and subset matrix
###############################################
outfile = file.path(outdir, "Voxel_Matrix.Rda")
load( file=outfile )
all.mat = mat
all.rs = rs
# ncut = 0.1
ncut = 10


template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)
dtemp = dim(temp)

t.t1 = file.path(tempdir, "sct1_unsmooth.nii.gz")
temp.t1 = readNIfTI(t.t1)


    
#############################################
# Overall desnity image
##############################################
noivh = demog$IVH_Dx_10 == 0 & demog$patientName != 131316
noivh.mat = mat[, which(noivh), drop=FALSE]

ivhvox = colSums(noivh.mat[ind, ])

rs = rowSums(noivh.mat)

all.ind = col.lists[["EVE_1"]]
all.ind = unlist(all.ind[ names(all.ind) != "Outside Brain Mask"])
V.noivh = sum(noivh.mat[all.ind,])
V = sum(all.mat[all.ind,])


csf.ind = col.lists[["EVE_1"]][["Background"]]
csfmask = array(0, dim = dtemp)
csfmask[csf.ind] = 1
csfmask = niftiarr(temp, csfmask)

sum(all.mat[csf.ind,])/V
sum(noivh.mat[csf.ind,])/V.noivh


csf.true = colSums(mat[csf.ind, ])
true = colSums(mat)
rat = csf.true/true

sum(noivh.mat[csf.ind, ])/sum(noivh.mat)
sum(all.mat[csf.ind, ])/sum(all.mat)

noivh.keep = which(rowSums(noivh.mat) > 0)
mean(noivh.mat[csf.ind,][, ])
#########################
# Of those with any hemorrhage, which percent are csf?
#########################
any.hem = which(rowSums(all.mat) > 10)
any.noivh.hem = which(rowSums(noivh.mat) > 6)

mean(csf.ind %in% any.hem)
mean(csf.ind %in% any.noivh.hem)



###########################
# Get mean prevalence for all non-zero CSF voxels
###########################
any.hem = rowSums(all.mat[csf.ind, ]) > 0
any.noivh.hem = rowSums(noivh.mat[csf.ind,]) > 0




mean(all.mat[csf.ind, ][any.hem,])
mean(noivh.mat[csf.ind, ][any.noivh.hem,])
mean(noivh.mat[csf.ind, ][any.hem,])

mean(all.mat[csf.ind, ])

mean(noivh.mat[csf.ind, ])
m = rowMeans(noivh.mat[csf.ind, ])
mean(m)
sd(m)
hist(colMeans(noivh.mat[csf.ind, ]))





noivh.csf.true = colSums(noivh.mat[csf.ind, ])
noivh.true = colSums(noivh.mat)
noivh.rat = noivh.csf.true/noivh.true

noivh.fac = factor(c("Has IVH", "No IVH")[noivh+1])
plot(rat ~ noivh.fac )


if (ncut < 1){
  ncut = floor(sum(noivh) * ncut)
}
noivh.mat = noivh.mat[rs > ncut, ]

runpct = rs / ncol(noivh.mat)

pngname = file.path(figdir, 
  paste0("Density_NoIVH.png"))

mbreak = ceiling(max(runpct) / 0.05) * 0.05
breaks = seq(0, mbreak, by=.05)
col.cut = alpha(div_gradient_pal(low="blue", 
                                 mid="red", 
                                 high="yellow")(
                                   seq(0, 1, 
                                    length=length(breaks)-1)
                                 ), .7)
pimg = niftiarr(temp, runpct)
clist = img_cut(pimg, breaks=breaks, include.lowest=FALSE)  
pimg[pimg == 0] = NA
cimg = clist$img
levs = clist$levs

plevs = levs
plevs = gsub("\\(", "", plevs)
plevs = gsub("\\]", "", plevs)
plevs = strsplit(plevs, ",")
plevs = lapply(plevs, as.numeric)
plevs = lapply(plevs, `*`, 100)
plevs = lapply(plevs, function(x){
  x[2] = x[2] - .01
  x
})
plevs = sapply(plevs, function(x){
  paste0(x[1], "-", x[2], "%")
})

png(pngname, res=600, height=7, width=7, units="in")
ortho2(temp.t1, pimg, col.y=col.cut, 
       ybreaks = breaks, 
       addlegend = TRUE,
       leg.x = 5, leg.y= 61, 
       legend=plevs, 
       leg.col=col.cut, leg.cex=1.5,
       leg.title = 
       "Percent of ICH Only Sample\n with Hemorrhage", 
       window = c(300, 1000))
dev.off()


pimg = niftiarr(temp, runpct)
dens.pct = area_pct(pimg, ind.list=lists[["EVE_1"]], 
                    keepall=TRUE)  
col.dens.pct = area_pct(pimg, ind.list= col.lists[["EVE_1"]], 
  keepall=TRUE)  

col.dens.tab = make.pvalimg(pimg, col.lists)[["EVE_1"]]
col.dens.tab = col.dens.tab[ order(col.dens.tab$nvox, 
                                   decreasing = TRUE), ]
col.dens.tab = col.dens.tab[ 1:10, ]

dens.tab = make.pvalimg(pimg, lists)[["EVE_1"]]
dens.tab = dens.tab[ order(dens.tab$nvox, 
                           decreasing = TRUE), ]
dens.tab = dens.tab[ 1:10, ]

save(dens.tab, dens.pct, 
     col.dens.pct, col.dens.tab, 
     file = file.path(figdir, 
                      "Table_NoIVH_Density_Results.Rda"))


rm(list=c("noivh", "noivh.mat", "ivhvox", "pimg"))


ind = col.lists[["EVE_1"]][["Outside Brain Mask"]]
ind.arr = array(0, dim=dtemp)
ind.arr[ind] = 1
ob = niftiarr(temp, ind.arr)

for (outcome in c("GCS", "NIHSS")){
  for (ncut in c(0.1, 10)){
    # outcome = "NIHSS"
    adder = paste0(outcome, "_")
    if (outcome == "NIHSS"){
      adder = ""
    }
    
    demog = all.demog
    if (outcome == "GCS") {
      demog$Y = demog$Enrollment_GCS_Add
    } else if (outcome == "NIHSS"){
      demog$Y = demog$Enrollment_NIHSS_Total
    } else {
      stop(paste0("Outcome ", outcome, " not implemented"))
    }
    
    
    xdemog = demog
    
    cc = complete.cases(demog$Y) & demog$IVH_Dx_10 == 0
    demog = demog[cc,]
    
    
    zform = ~ Age + Sex + Diagnostic_ICH
    # Z = model.matrix(object = zform, data = demog)
    # Z = model.matrix(object = zform, data = demog)
    
    ###############################################
    # Subset matrix
    ###############################################  
    mat = all.mat
    rs = all.rs
      
    mat = mat[, which(cc), drop=FALSE]
    pct = FALSE
    if (ncut < 1){
      pct = TRUE
      ncut = floor(nrow(demog) * ncut)
      rs = rowSums(mat)
    }
    dim(mat)
    sum(rs > ncut)
    mat = mat[rs > ncut, ]
    
    for (icut in 0:3){
      pngname = file.path(figdir, 
        paste0("Density_NoIVH_Above_", icut, 
         "_Patients_", outcome, ".png"))
      if (!file.exists(pngname) & pct){
        img = array(0, dim= dtemp)
        img[ rs > ncut] = 1
        img = niftiarr(temp, img)
        #   mask.overlay(temp, img)
        img[ !is.na(img) ] = rs / ncol(mat)
        #   img[ rs == 0 ] = NA
        img[ rs <= icut ] = NA
        img = cal_img(img)
        
        png(pngname)
        ortho2(temp, img, col.y= hotmetal()[10:64])
        dev.off()
        
      }
    
    }


    
    if (pct){
      runpct = rs / ncol(mat)
      mbreak = ceiling(max(runpct) / 0.05) * 0.05
      breaks = seq(0, mbreak, by=.05)
      col.cut = alpha(div_gradient_pal(low="blue", 
                                       mid="red", 
                                       high="yellow")(
                                         seq(0, 1, 
                                          length=length(breaks)-1)
                                       ), .7)
      pimg = niftiarr(temp, runpct)
      dens.pct = area_pct(pimg, ind.list=lists[["EVE_1"]], 
                          keepall=TRUE)  
      col.dens.pct = area_pct(pimg, ind.list= col.lists[["EVE_1"]], 
                              keepall=TRUE)  
      
      col.dens.tab = make.pvalimg(pimg, col.lists)[["EVE_1"]]
      col.dens.tab = col.dens.tab[ order(col.dens.tab$nvox, 
                                       decreasing = TRUE), ]
      col.dens.tab = col.dens.tab[ 1:10, ]
          
      dens.tab = make.pvalimg(pimg, lists)[["EVE_1"]]
      dens.tab = dens.tab[ order(dens.tab$nvox, 
                                         decreasing = TRUE), ]
      dens.tab = dens.tab[ 1:10, ]
      
      save(dens.tab, dens.pct, 
           col.dens.pct, col.dens.tab, 
           file = file.path(figdir, 
                            paste0(outcome, "_NoIVH_Density_keep",
                                   ncut,"_Results.Rda")
           )
      )
      
      clist = img_cut(pimg, breaks=breaks, include.lowest=FALSE)  
      pimg[pimg == 0] = NA
      cimg = clist$img
      levs = clist$levs
      
      plevs = levs
      plevs = gsub("\\(", "", plevs)
      plevs = gsub("\\]", "", plevs)
      plevs = strsplit(plevs, ",")
      plevs = lapply(plevs, as.numeric)
      plevs = lapply(plevs, `*`, 100)
      plevs = lapply(plevs, function(x){
        x[2] = x[2] - .01
        x
      })
      plevs = sapply(plevs, function(x){
        paste0(x[1], "-", x[2], "%")
      })
      

      
      pngname = file.path(figdir, paste0("Density_Image_NoIVH_", 
        outcome, ".png"))
      png(pngname, res=600, height=7, width=7, units="in")
      ortho2(temp, pimg, col.y=col.cut, 
             ybreaks = breaks, 
             addlegend = TRUE,
             leg.x = 5, leg.y= 60, 
             legend=plevs, 
             leg.col=col.cut, leg.cex=1.5,
             leg.title = 
             "Percent of ICH Only Sample\n with Hemorrhage")
      dev.off()
      pimg[ pimg < 0.05] = NA

      pngname = file.path(figdir, paste0("Density_Image_NoIVH_", 
        outcome, "_Over5pct.png"))
      png(pngname, res=600, height=7, width=7, units="in")
      ortho2(temp, pimg, col.y=col.cut, 
             ybreaks = breaks, 
             addlegend = TRUE,
             leg.x = 5, leg.y= 60, 
             legend=plevs, 
             leg.col=col.cut, leg.cex=1.5,
             leg.title = 
             "Percent of ICH Only Sample\n with Hemorrhage")
      dev.off()    
    }
    ###############################################
    # Cross validation folds
    ###############################################
    
    runmat = t(mat)
    class(runmat) = "numeric"
    
    ###############################################
    # Get Y values
    ###############################################  
    runY = demog$Y
    print(outcome)
    mytime = system.time({
      mods = fast_lm(runY, X=runmat, Z = NULL, 
                     spot.check = TRUE, ncheck = 100, 
                     verbose=FALSE)
    })  
    
    pvals = lm(runmat ~ runY)
    mods$AIC = extractAIC(mods)
    
    icut = 0.05
    for (icut in c(0.01, 0.001, 0.05, 1000, 2000, 3000)){
      img = array(NA, dim= dtemp)
      
      if (icut < 1){
        ind = mods$p.val <= icut
      } else {
        ind  = rank(mods$p.val) <= icut        
      }
      
      #       if (outcome == "NIHSS"){
      #         pval = 0.01
      #         ind = mods$p.val <= pval
      #       }
      #       if (outcome == "GCS"){
      #         ord = 1000
      #         ind  = rank(mods$p.val) <= ord
      #       }
      img[ rs > ncut ][ind] = 1
      img[is.na(img)] = 0
      img = niftiarr(temp, img)
      img = cal_img(img)
      
      pngname = file.path(figdir, paste0("ROI_", icut, "_keep", 
        ncut, "_", outcome, ".png"))
      png(pngname)
      mask.overlay(temp, img)  
      dev.off()
      
      symm.img = lr_symm(img)
      pngname = file.path(figdir, paste0("Symm_ROI_", icut,  
        "_keep", ncut, "_", outcome, ".png"))
      png(pngname)
      mask.overlay(temp, symm.img)  
      dev.off()  
      
      symm.ind = which(symm.img > 0)
      
      rr = rowMeans(runmat[, ind, drop=FALSE])
      demog[, paste0("pct")] = rr
      
      all.ind = which(img > 0)
      all.rr = colMeans(all.mat[all.ind, , drop=FALSE])
      xdemog[, paste0("pct")] = all.rr
      
      all.rr = colMeans(all.mat[symm.ind, , drop=FALSE])
      xdemog[, paste0("symm_pct")] = all.rr
      
      sub = xdemog[ all.rr == 0, ]
      
      pvalimg.pct = area_pct(img, ind.list=lists[["EVE_1"]], 
        keepall=TRUE)  
      pvalimg.tab = make.pvalimg(img, lists)[["EVE_1"]]
      pvalimg.tab = pvalimg.tab[ order(pvalimg.tab$nvox, 
        decreasing = TRUE), ]
      pvalimg.tab = pvalimg.tab[ 1:10, ]
      
      col.pvalimg.pct = area_pct(img, ind.list=
        col.lists[["EVE_1"]], keepall=TRUE)  
      
      col.pvalimg.tab = make.pvalimg(img, col.lists)[["EVE_1"]]
      col.pvalimg.tab = col.pvalimg.tab[ 
        order(col.pvalimg.tab$nvox, decreasing = TRUE), ]

      col.pvalimg.tab = col.pvalimg.tab[ 1:10, ]
      
      pvalimg.tab_symm = 
        make.pvalimg(symm.img, lists)[["EVE_1"]]
      col.pvalimg.tab_symm = 
        make.pvalimg(symm.img, col.lists)[["EVE_1"]]
      
      save(pvalimg.pct, pvalimg.tab, 
           col.pvalimg.pct, col.pvalimg.tab, 
           pvalimg.tab_symm, col.pvalimg.tab_symm, 
           img,
           file = file.path(figdir, 
            paste0(outcome, "_NoIVH_", icut,  
             "_keep", ncut,"_Results.Rda")
                            )
           )
    }
  }
} 




