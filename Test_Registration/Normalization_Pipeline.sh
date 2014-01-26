#!/bin/bash 
cd /dexter/disk2/smart/stroke_ct/ident/programs/Test_Registration

qsub -cwd -l mem_free=5G,h_vmem=7G -N MakeNifti make_Registration_nifti_array.sh

qsub -cwd -hold_jid MakeNifti -N Skull_Strip Test_Registration_Run_Skull_Strip_Array.sh

# rootdir="/dexter/disk2/smart/stroke_ct/ident"
# basedir="${rootdir}/Test_Registration"
# refdir="${basedir}/RawNIfTI/Skull_Stripped/"
# outdir="${basedir}/FLIRT/Skull_Stripped"
# mkdir -p "$outdir"

# cp $refdir/*.nii* $outdir/

# qsub -cwd -hold_jid MakeNifti -l mem_free=5G,h_vmem=7G -N reorient reorient_files.sh

qsub -cwd -hold_jid MakeNifti -l mem_free=5G,h_vmem=7G -N ClinNorm \
	run_clinical_normalization.sh