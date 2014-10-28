#!/bin/bash 

cd $dex/programs/Segmentation; 

qcwd -l mem_free=12G,h_vmem=13G -N REG CT_Register_to_Template_SyN.sh

qcwd -l mem_free=30G,h_vmem=33G -N INDIV -hold_jid_ad REG \
Run_Indiv_Model.sh

Rbatch Get_True_Vol.R -N TRUTH -hold_jid INDIV \
-l mem_free=10G,h_vmem=11G


### Maxvmem = 31G
qcwd -l mem_free=30G,h_vmem=31G -N AGGDATA \
-hold_jid INDIV Create_Aggregate_Data.sh

qcwd -l mem_free=10G,h_vmem=12G -N INCLUDE \
-hold_jid INDIV -hold_jid_ad AGGDATA Inclusion_Criteria_Results.sh

qcwd -l mem_free=15G,h_vmem=17G -N PLOT \
-hold_jid AGGDATA Plot_Agg_Relationship.sh


qcwd -l mem_free=30G,h_vmem=31G -N AGG -hold_jid INDIV \ 
-hold_jid AGGDATA Aggregate_Model.sh

qcwd -l mem_free=10G,h_vmem=12G -N COLL -hold_jid INDIV \
-hold_jid AGG Collapse_Predictor_Models.sh

qcwd -l mem_free=30G,h_vmem=31G -N CUT \
-hold_jid COLL Get_Model_Cutoffs.sh

qcwd -l mem_free=20G,h_vmem=21G -N SCUT \
-hold_jid COLL Get_Smooth_Model_Cutoffs.sh

qcwd -l mem_free=20G,h_vmem=21G -N PRED \
-hold_jid CUT Predict_From_Model.sh

Rbatch Get_Results.R -l mem_free=8G,h_vmem=10G \
-N RES -hold_jid PRED


