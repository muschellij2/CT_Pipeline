cd /dexter/disk2/smart/stroke_ct/ident/programs/Test_Registration

qsub -cwd -l mem_free=5G,h_vmem=7G -N MakeNifti make_Registration_nifti_array.sh

qsub -cwd -hold_jid MakeNifti -l mem_free=5G,h_vmem=7G -N ClinNorm run_clinical_normalization.sh