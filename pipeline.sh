#!/bin/bash 

rootdir="/dexter/disk2/smart/stroke_ct/ident"
basedir="${rootdir}/Test_5"
progdir="${rootdir}/programs"
nfolds=75;
#"\\d\\d\\d-(\\d|)\\d\\d\\d"
nid=`find $basedir -mindepth 1 -maxdepth 1 -type d -regex '.*[0-9][0-9][0-9]-[0-9]?[0-9][0-9][0-9]$' | wc -l`
# find $basedir -regex '{[0-9]}'

cd $progdir

### Make NIfTI Files
qsub -cwd -N nifti_maker -l mem_free=2G,h_vmem=8G \
make_nifti_array.sh

### Brain extraction and mask generation
qsub -cwd -hold_jid nifti_maker -N Skull_Strip \
Run_Skull_Strip_Array.sh

### create individual data sets
qsub -cwd -hold_jid nifti_maker \
-hold_jid_ad Skull_Strip -N make_df create_id_df.sh

### create collapsed data set
qsub -cwd -hold_jid make_df -N coll_data collapse_data.sh

### create models
qsub -cwd -hold_jid Skull_Strip \
-hold_jid_ad make_df -N predict_clot \
-l mem_free=30G,h_vmem=20G predict_clot_array.sh

### create collapsed model
qsub -cwd -hold_jid predict_clot -hold_jid coll_data \
-l mem_free=30G,h_vmem=25G -N coll_mod collapse_models.sh

### create collapsed model
qsub -cwd -l mem_free=20G,h_vmem=40G \
-hold_jid coll_mod -N cmod_roc collapse_models_ROC.sh

### create collapsed model
qsub -cwd -hold_jid coll_mod -N train_data First_Train_Data.sh

### create collapsed model
qsub -cwd -hold_jid train_data -N predmods \
-l mem_free=1G,h_vmem=2G multiple_prediction_models.sh

### create collapsed model
qsub -cwd -hold_jid predmods -N aggmod \
-l mem_free=1G,h_vmem=2G multiple_predictions.sh


### create collapsed model
qsub -cwd -hold_jid aggmod -N auc \
-l mem_free=12G,h_vmem=13G multiple_AUC.sh







