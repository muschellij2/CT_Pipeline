#################################
# Sub_analyses of Groups
# Author: John Muschelli
#################################
rm(list=ls())
library(cttools)
library(fslr)
library(oro.nifti)
library(scales)
library(RColorBrewer)
library(data.table)
library(ggplot2)
library(grid)
library(plyr)
homedir = "/Applications"
rootdir = "~/CT_Registration"
basedir = file.path(rootdir, "data")
outdir = basedir
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
  basedir = file.path(rootdir, "Registration")
  outdir = file.path(basedir, "results")
}
progdir = file.path(rootdir, "programs")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")


template = file.path(tempdir, "scct_unsmooth.nii.gz")
temp = readNIfTI(template)

t.t1 = file.path(tempdir, "sct1_unsmooth.nii.gz")
temp.t1 = readNIfTI(t.t1)

shell.img = temp.t1
shell.img [is.na(shell.img ) | !is.na(shell.img )] = 0 
shell.img = cal_img(shell.img)
shell.img = zero_trans(shell.img)

whichdir = "reoriented"
outcome = "NIHSS"
adder = paste0(outcome, "_")
if (outcome == "NIHSS"){
  adder = ""
}

get.id = function(x){
  ss = strsplit(x, "_")
  ss = sapply(ss, head, 1)
  ss = gsub(".*(\\d\\d\\d-.*)", "\\1", ss)
  ss
}

id_to_pname = function(x){
  as.numeric(gsub("-", "", x))
}

demog = read.csv(file=file.path(basedir, "Demog_NIHSS_Mask.csv"), 
 stringsAsFactors=FALSE)
demog$Base_ICH_10 = demog$Diagnostic_ICH /10

if (outcome == "GCS") {
  demog$Y = demog$Enrollment_GCS_Add
} else if (outcome == "NIHSS"){
  demog$Y = demog$Enrollment_NIHSS_Total
} else {
  stop(paste0("Outcome ", outcome, " not implemented"))
}
 

demog$Clot_Location_RC = gsub("Palidus", "Pallidus", 
  demog$Clot_Location_RC )
demog$Clot_Location_RC = factor(demog$Clot_Location_RC, 
  levels= c("Lobar", "Globus Pallidus", "Putamen", "Thalamus"))
demog$LOC = demog$Clot_Location_RC

demog$GCS = ifelse(demog$Enrollment_GCS_Add <= 12, "3-12" ,"13-15")
gcs.levs = c("3-12" ,"13-15")
demog$GCS = factor(demog$GCS, levels=gcs.levs)

demog$Location = factor(demog$Lobar_BG_vs, 
  levels=c("Deep", "Lobar"))

# Z = model.matrix(object = zform, data = demog)
# Z = model.matrix(object = zform, data = demog)

outfile = file.path(outdir, "Voxel_Matrix.Rda")
load( file=outfile )

keep = rs > 0
check.ids = id_to_pname(get.id(colnames(mat)))
stopifnot(all(check.ids == demog$patientName))

out = "GCS"
# out = "Location"

cc = complete.cases(demog[, out])

out.levs = levels(demog[, out])

cc.demog = demog[cc, ]

submat = mat[keep, cc]

ind = sapply(out.levs, function(x){
  which(cc.demog[, out] == x)
})

props = sapply(ind, function(x){
  rowMeans(submat[, x])
})

make.img = function(x) {
  img = shell.img
  img@.Data[rs > 0] = x 
  img = cal_img(img)
  return(img)
}


img_cut = function(img, breaks, ...){
  cuts = cut(img, breaks=breaks, ...)
  # cuts = factor(cuts, levels)
  levs = levels(cuts)
  cuts = as.numeric(cuts)
  # res.p[ rs > ncut ] = cuts
  img@.Data = array(cuts, dim=dim(img))
  return(list(img=img, levs=levs))
}


pimgs = alply(props, 2, make.img)
names(pimgs) = colnames(props)

ranges = sapply(pimgs, range)
ranges = c(min(ranges), max(ranges))

breaks = seq(ranges[1], ranges[2], by=.05)

col.cut = alpha(div_gradient_pal(low="blue", 
  mid="red", 
  high="yellow")(
  seq(0, 1, length=length(breaks)-1)
  ), .7)

img = check_nifti(temp.t1)
window = c(300, 1000)
img@cal_min = window[1]
img@cal_max = window[2]
img[img < window[1]] = window[1]
img[img >= window[2]] = window[2]

i = 1

for(i in seq_along(out.levs)){

  pimg = pimgs[[i]]

  img.mask = check_nifti(pimg)

  img.mask[img.mask <= 0] = NA

  clist = img_cut(pimg, breaks=breaks, include.lowest=FALSE)
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


  out.lev = out.levs[i]
  fname = file.path(outdir, 
    paste0(out, "_", out.lev, "_3D_Histogram.png"))
  png(fname, res=600, height=7, width=7, units="in", type="cairo")
  ortho2(img, img.mask, col.y=col.cut, 
    ybreaks = breaks, 
    addlegend = TRUE,
    leg.x = 13, leg.y= 60, 
    legend=levs, 
    leg.col=col.cut, leg.cex=1.5,
    leg.title = paste0("Proportion of ", out, ": ", out.lev, 
      "\n with Hemorrhage"))
  dev.off()
  print(fname)
}