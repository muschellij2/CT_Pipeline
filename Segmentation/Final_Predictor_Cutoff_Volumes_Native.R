###########################################################
## This calculates the volumes for the Rigid back trans
## formed data
##
## Author: John Muschelli
## Last updated: May 20, 2014
###########################################################
###########################################################
rm(list=ls())
library(plyr)
library(cttools)
library(fslr)
library(ROCR)
library(matrixStats)
library(extrantsr)
library(ANTsR)
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
segdir = file.path(progdir, "Segmentation")
source(file.path(segdir, "performance_functions.R"))


template.file = file.path(tempdir, "scct_unsmooth.nii.gz")
ss.tempfile = file.path(tempdir, "Skull_Stripped",
    "scct_unsmooth_SS_First_Pass_0.1.nii.gz")

#### load voxel data
outfile = file.path(outdir, "Voxel_Info.Rda")
load(file=outfile )

outfile = file.path(outdir, 
    "111_Filenames_with_volumes_stats.Rda")
load(file = outfile)

regs = "Rigid"
scen = expand.grid(ttype = regs,
    interpolator = c("Linear"),
    stringsAsFactors = FALSE)
iscen = 1
# scen = scen[2,, drop=FALSE]
# addons = paste0(scen$ttype, "_", scen$interpolator) 

# fdf[, addons] = NA

# types = c("_zval2", '_zval2_medztemp')
types = "_zval2"
# , "_zval2"
# "_include_all", 
type = types[1]

# for (iscen in seq(nrow(scen))){
# addon = addons[iscen]

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
correct = gsub("^_", "", img_ext)
adder = switch(correct, 
    "none"= "",
    "N3"="_N3",
    "N4" = "_N4",
    "N3_SS" = "_N3_SS",
    "N4_SS" = "_N4_SS", 
    "SyN" = "_SyN",
    "SyN_sinc" = "_SyN_sinc",
    "Rigid" = "_Rigid",
    "Affine" = "_Affine",
    "Rigid_sinc" = "_Rigid_sinc",
    "Affine_sinc" = "_Affine_sinc")    


mod.filename = file.path(outdir, 
    paste0("Collapsed_Models", adder, ".Rda"))
load(mod.filename)
cns = colnames(res)

    
cut.filename = file.path(outdir, 
paste0("Model_Cutoffs", adder, ".Rda"))

load(file=cut.filename)

scut.filename = file.path(outdir, 
paste0("Smooth_Model_Cutoffs", adder, ".Rda"))

load(file=scut.filename)
names(all.sdice.cuts) = names(all.scuts) = cns


makedir = sapply( fdf$outdir, function(x) {
    if (!file.exists(x)){
        dir.create(x, showWarnings =FALSE)
    }
})
fdf$ss = gsub("_Mask", "", fdf$mask)

fdf$hdr = file.path(fdf$iddir, "Sorted",
    paste0(nii.stub(fdf$img, bn=TRUE), 
        "_Header_Info.Rda"))

iimg <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(iimg)) iimg = 71

cns = c("mod_agg", "gam", "rf")

cutter = c("", "_dice")
icut = cutter[1]


vals = c(outer(types, cutter, paste0))
vals = c(outer(vals, c("_pred"), paste0))
vals = c(outer(cns, vals, paste0))

vals = paste0("vol_", vals)
xvals = vals

vals = paste0("vol_", c(outer(cns, types, paste0)), "_raw")
vals = c(vals, xvals)

xvals = paste0("vol_", c(outer(cns, types, paste0)), "_1")
vals = c(vals, xvals)

fdf[, vals] = NA

reses = matrix(NA, nrow=nrow(fdf), 
    ncol=prod(length(cns), length(types), length(cutter)))
colnames(reses) = c(outer(
    outer(cns, types, paste0), 
    cutter, 
    paste0)
)

dice.mod.saccs = reses
dice.mod.ssens = reses
dice.mod.sspec = reses
dice.mod.sdice = reses
rm(list="reses")

# ss = c(outer( paste0(cns, adder), types, paste0))
# ss_smooth = c(outer( paste0(cns, adder, 
    #"_smoothed", adder), 
#     types, paste0))


for (iimg in seq(nrow(fdf))){

    x = fdf[iimg,]    
    vol = x$truevol

    ####################################
    ## Run both with the Skull Stripped and 
    ## not skull stripped
    ####################################
    x$preddir = x$outdir
    x$outdir = file.path(x$iddir, outputdir)
    hdrload = load(x$hdr)

    ########
    #!!! need to load in ROI
    ########    
    roi = readNIfTI(x$roi, reorient=FALSE) > 0.5


    # cn = "gam"
    cn = cns[1]
    for (cn in cns){
        for (type in types){

        # ss = c(outer( paste0(cn, adder), type, paste0))
            ss_smooth = c(outer( paste0(cn, adder, 
                "_smoothed", adder), 
                type, paste0))

            outimg = nii.stub(x$img, bn=TRUE)
            outimg = file.path(x$preddir, 
                paste0(outimg, "_", ss_smooth))

            outfile = paste0(outimg, "_native.nii.gz")

            native.roi = cal_img(readNIfTI(outfile, 
                    reorient=FALSE))
            vol.raw = get_roi_vol(native.roi, 
                dcmtables)$truevol

            val = paste0("vol_", cn, type, "_raw")
            fdf[iimg, val] = vol.raw

            vol.1 = get_roi_vol(native.roi >= 0.9, 
                dcmtables)$truevol

            val = paste0("vol_", cn, type, "_1")
            fdf[iimg, val] = vol.1

            print(outfile)
            cat("True Volume is \n")
            print(vol)
            cat("Raw Volume is \n")
            print(vol.raw)

            for (icut in cutter){
                runcol = paste0(cn, type, icut)

                if (icut == ""){
                    cutval = all.scuts[cn]
                }
                if (icut == "_dice"){
                    cutval = all.sdice.cuts[cn]
                }   

                native.pred = cal_img(native.roi > cutval)
                
                tab = extrantsr:::my.tab(native.pred, 
                    roi)

                dice.mod.saccs[iimg, runcol] = 
                    get.acc(tab)
                dice.mod.ssens[iimg, runcol] = 
                    get.sens(tab)
                dice.mod.sspec[iimg, runcol] = 
                    get.spec(tab)
                dice.mod.sdice[iimg, runcol] = 
                    get.dice(tab)

                vol.pred = get_roi_vol(native.pred, 
                    dcmtables)$truevol
                cat(paste0("Cut Volume for ", 
                    icut, " cutoff is\n"))
                print(vol.pred)
                val = paste0("vol_", cn, 
                    type, icut, "_pred")
                fdf[iimg, val] = vol.pred
            }

        }
    }
    # }
    # 
    print(iimg)
}

outfile = file.path(outdir, 
    paste0("Cutoff_predictions.Rda"))
save(fdf, 
    dice.mod.saccs, 
    dice.mod.ssens, 
    dice.mod.sspec, 
    dice.mod.sdice,
    file = outfile)

    # print(iscen)
# }

