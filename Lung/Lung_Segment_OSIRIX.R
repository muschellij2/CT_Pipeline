#################################
# Lung Segmentation for OsiriX data
#################################
rm(list=ls())
library(cttools)
library(fslr)
library(plyr)
library(R.utils)
options(matlab.path='/Applications/MATLAB_R2013b.app/bin')

rootdir = "~/Lung/data"
rootdir = path.expand(rootdir)

iddirs = list.dirs(path=rootdir, full.names=TRUE, recursive=FALSE)
ids = basename(iddirs)

iid <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if (is.na(iid)) iid <- 5

# for (iid in seq_along(ids)){
    basedir = iddirs[iid]

    niis = list.files(path=basedir, 
        pattern= ".nii.gz", 
        full.names=TRUE, recursive=FALSE)
    outdir = file.path(basedir, "Lung")
    dir.create(outdir, showWarnings=FALSE)
    outstubs = file.path(outdir, 
        nii.stub(niis, bn=TRUE))
    outfiles = paste0(outstubs, "_Lung_Mask")
    lfiles = paste0(outstubs, "_Lung")
    hfiles = paste0(outstubs, "_Human")
    # exfiles = paste0(outfiles, "_Regioned")
    # smfiles = paste0(outstubs, "_Smooth")
    ifile = 1
    ### do lung segmentation
    for (ifile in seq_along(niis)){
        ofile = outfiles[ifile]
        hfile = hfiles[ifile]
        lfile = lfiles[ifile]
        img = niis[ifile]
        img = check_nifti(img)
        # exfile = exfiles[ifile]
        himg = CT_human_mask(img, retimg=TRUE, 
            outfile = paste0(hfile, ".nii.gz"))
        lung_mask = CT_lung_mask(himg, human = FALSE,
            outfile = paste0(ofile, ".nii.gz"), topN =1)

        lung = img
        lung[lung_mask == 0] = -1024
        lung = cal_img(lung)
        writeNIfTI(lung, filename = lfile)


        # mask = segment_lung(img, outfile = ofile,
        #     smooth.outfile = smfile)
        # spm_bwlabel(infile=paste0(ofile, ".nii.gz"), 
        #     outfile = paste0(exfile, ".nii.gz"), topN = 2)
    }
    print(iid)
    