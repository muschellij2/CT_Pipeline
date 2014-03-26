#!/bin/bash 

rootdir="/dexter/disk2/smart/stroke_ct/ident"
basedir="${rootdir}/Registration"
progdir="${rootdir}/programs"
#"\\d\\d\\d-(\\d|)\\d\\d\\d"
nid=`find $basedir -mindepth 1 -maxdepth 1 -type d -regex '.*[0-9][0-9][0-9]-[0-9]?[0-9][0-9][0-9]$' | wc -l`
# find $basedir -regex '{[0-9]}'

#### need to unzip the ROIs then run Unzip_ROIs.R
#### then run Registration_Pipeline wih ROI_data as study and ROIformat
##### = TRUE

#NIfTI Files

cd $progdir/ROI

### Copy over files and turn to .nii for SPM
qsub -cwd -N Copying -l mem_free=2G,h_vmem=8G ROI_Copying_Reorienting.sh

cd $progdir

### Spatially Normalize
qsub -cwd -hold_jid Copying -N NORMALIZE -l mem_free=20G,h_vmem=22G \
Reorient_Normalize_Subsample_Overlap.sh

### Get overlap of atlases
qsub -cwd -hold_jid NORMALIZE -N OVERLAP -l mem_free=20G,h_vmem=22G \
Atlas_Overlap_Report.sh 

### create REport
qsub -cwd -hold_jid OVERLAP -N REPORT -l mem_free=12G,h_vmem=14G \
100_Reg_Report.sh

###
cd $progdir/ROI

qsub -cwd -l mem_free=20G,h_vmem=21G -hold_jid NORMALIZE \
 -N VOXMAT ROI_Voxel_Matrix.sh

qsub -cwd -l mem_free=20G,h_vmem=21G -N REGRESS \
-hold_jid VOXMAT NIHSS_Regress.sh







