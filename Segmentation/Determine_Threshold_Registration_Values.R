####################################################################
## This code is for Thresholding Registered data
##
## Author: John Muschelli
## Last updated: May 20, 2014
####################################################################
####################################################################
rm(list=ls())
library(plyr)
library(ggplot2)
library(matrixStats)
library(reshape2)
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

#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )

outfile = file.path(outdir, "111_Filenames_with_volumes.Rda")
load(file = outfile)

# regs = c("Affine", "Rigid", "SyN")
regs = "Rigid"
scen = expand.grid(ttype = regs,
    interpolator = c("Linear", "LanczosWindowedSinc"),
    stringsAsFactors = FALSE)
iscen = 1
addons = paste0(scen$ttype, "_", scen$interpolator) 
addons.diff = paste0(addons, ".diff")
addons.adiff = paste0(addons, ".adiff")
addons.est = paste0(addons, ".est")

fdf[, addons.adiff] = fdf[, addons.diff] = fdf[, addons] = NA
fdf[, addons.est] = fdf[, addons.adiff] 

N = nrow(scen)
cuts = seq(0.01, 1, by=0.001)
mat = array(NA, dim=c(length(cuts), nrow(fdf), N))


for (iscen in seq(N)){
    
    addon = addons[iscen]
    addon.diff = addons.diff[iscen]
    addon.adiff = addons.adiff[iscen]
    addon.est = addons.est[iscen]
    
    # ttype = "Rigid"
    ttype = scen$ttype[iscen]
    # ttype = "SyN"
    # interpolator = c("Linear", "LanczosWindowedSinc")
    interpolator = scen$interpolator[iscen]
    int_ext = switch(interpolator,
        "Linear" = "",
        "LanczosWindowedSinc" = "_sinc")
    outputdir = paste0(ttype, "_Registered")
    img_ext = paste0("_", ttype, int_ext)


    # template.file = file.path(tempdir, "scct_unsmooth.nii.gz")
    # ss.tempfile = file.path(tempdir, "Skull_Stripped",
    #     "scct_unsmooth_SS_First_Pass_0.1.nii.gz")


    makedir = sapply( fdf$outdir, function(x) {
    	if (!file.exists(x)){
    		dir.create(x, showWarnings =FALSE)
    	}
    })
    fdf$ss = gsub("_Mask", "", fdf$mask)
    iimg = 21

    for (iimg in seq(nrow(fdf))){
        vol = fdf$truevol[iimg]

        x = fdf[iimg,]

        stubfile = function(x, d = NULL, ext = ""){
          b = nii.stub(x, bn=TRUE)
          b = paste0(b, ext)
          file.path(d, b)
        }

        ofile = stubfile(x$ss, d = x$outdir, ext = img_ext)
        outprefix = stubfile(x$ss, d = x$outdir)

        roi.ofile = paste0(ofile, "_ROI.nii.gz")
        mask.ofile = paste0(ofile, "_Mask.nii.gz")
        ofile = paste0(ofile, ".nii.gz")
        print(ofile)
        binary = c(roi.ofile, mask.ofile)
        files = c(ofile, binary)
        ex = all(file.exists(files))
        # if (!ex){
        roi =readNIfTI(roi.ofile, 
            reorient=FALSE)
        vals = roi[roi > 0]
        cuts = seq(0.01, 1, by=0.001)
        vres = prod(pixdim(roi)[2:4])
        stopifnot(vres == 1)
        est = sapply(cuts, function(x) sum(vals >= x)/1000)
        vdiff = est - vol
        adiff = abs(est - vol)
        best.est = which.min(adiff)
        best.cutoff = cuts[best.est]
        fdf[iimg, addon] = best.cutoff
        outfile = file.path(x$iddir, "Predictors", 
            paste0("Volume_cutoff_", addon, ".Rda"))
        save(best.cutoff, vdiff, adiff, cuts, 
            est, best.est, file = outfile)        

        ####################################
        ## Run both with the Skull Stripped and not skull stripped
        ####################################
        xcut = load(file = outfile) 

        mat[, iimg, iscen] = vdiff
        fdf[iimg, addon.est] = est[best.est]
        fdf[iimg, addon] = best.cutoff
        fdf[iimg, addon.diff] = vdiff[best.est]
        fdf[iimg, addon.adiff] = adiff[best.est]
        # fdf[iimg, addon.diff] = vdiff[best.est]
        print(iimg)
    }
    print(iscen)
}

lmat = vector(mode="list", length=N)
for (i in seq(N)) lmat[[i]] = mat[,,i]
lmat = lapply(lmat, function(x){
    colnames(x) = fdf$img
    m = melt(x)
    colnames(m) = c("row", "img", "value")
    m = cbind(m, cuts)
})

lmat = mapply(function(x, cn){
    x$scen = cn
    x
}, lmat, addons, SIMPLIFY=FALSE)

lmat = do.call("rbind", lmat)
lmat$row = NULL




cmed = colMedians(as.matrix(fdf[, addons]))
cm = colMeans(as.matrix(fdf[, addons]))
colMins(as.matrix(fdf[, addons]))

(cmed.diff = colMedians(as.matrix(fdf[, addons.diff])))
(cm.diff = colMeans(as.matrix(fdf[, addons.diff])))
(cmax.diff = colMaxs(as.matrix(fdf[, addons.diff])))

(cmed.adiff = colMedians(as.matrix(fdf[, addons.adiff])))
(cm.adiff = colMeans(as.matrix(fdf[, addons.adiff])))
(cmax.adiff = colMaxs(as.matrix(fdf[, addons.adiff])))

gli = ggplot(lmat, aes(x=cuts, y=value, colour=img)) + 
    guides(colour=FALSE) + geom_line() + facet_wrap(~ scen)

cut0.5 = lmat[ lmat$cuts == 0.5,]

cut0.5 = merge(cut0.5, fdf[, c("img", "varslice")], 
    by="img")

cms = round(cm, 3)
cms = data.frame(cuts = cms, scen = names(cms),
    stringsAsFactors = FALSE)
cutter = merge(lmat, cms, all.y=TRUE)

cutter = merge(cutter, fdf[, c("img", "varslice")], 
    by="img")

pdfname = file.path(outdir, "Cutoff_Histograms.pdf")
pdf(pdfname)
     
    g = ggplot(fdf, aes(x= truevol, 
        colour = factor(varslice))) + 
        geom_point() + geom_smooth(se=FALSE) + 
        geom_abline(xintercept = 0, 
        slope = 1) + guides(colour = FALSE)
    g + aes(y = Affine_LanczosWindowedSinc.est)
    g + aes(y = Affine_Linear.est)
    g + aes(y = SyN_LanczosWindowedSinc.est)
    g + aes(y = SyN_Linear.est)
    g + aes(y = Rigid_Linear.est)
    g + aes(y = Rigid_LanczosWindowedSinc.est)


    print(gli)

    g = ggplot(fdf, aes(fill = factor(varslice))) + 
    geom_histogram(position="identity", alpha = 0.5) + 
    guides(fill = FALSE) + xlim(c(0, 1))
    for (iadd in addons) {
        print({g + aes_string(x = iadd) + ggtitle(iadd)})
    }
    h05 = ggplot(cut0.5, aes(x=value, fill= factor(varslice))) + 
        guides(fill=FALSE) + 
        geom_histogram(position="identity", alpha = 0.5) + 
        facet_wrap(~ scen)
    print(h05 )
    print({h05 %+% 
        cut0.5[ grep("Rigid", cut0.5$scen), ]
    })
    print({h05 %+% cutter })    

dev.off()
