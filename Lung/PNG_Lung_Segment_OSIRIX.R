#################################
# Lung Segmentation for OsiriX data
#################################
rm(list=ls())
library(cttools)
library(fslr)
library(plyr)
options(matlab.path='/Applications/MATLAB_R2016a.app/bin')


rootdir = "~/Lung/data"
rootdir = path.expand(rootdir)

iddirs = list.dirs(path=rootdir, full.names=TRUE, recursive=FALSE)
ids = basename(iddirs)

iid <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if (is.na(iid)) iid <- 3

# for (iid in seq_along(ids)){
    basedir = iddirs[iid]

    niis = list.files(path=basedir, 
        pattern= ".nii.gz", 
        full.names=TRUE, recursive=FALSE)
    outdir = file.path(basedir, "Lung")
    dir.create(outdir, showWarnings=FALSE)
    outstubs = file.path(outdir, 
        nii.stub(niis, bn=TRUE))
    outfiles = paste0(outstubs, "_Lung")
    smfiles = paste0(outstubs, "_Smooth")
    ifile = 1
    ### do lung segmentation
    for (ifile in seq_along(niis)){
        ofile = outfiles[ifile]
        img = check_nifti(niis[ifile])
        opng = paste0(ofile, ".png")
        ofile = check_nifti(ofile)
        xyz= ceiling(cog(ofile))
        png(opng, type="cairo")
            mask.overlay(img, ofile, window=c(0, 200), xyz=xyz)
        dev.off()
    }    
    print(iid)

# }
