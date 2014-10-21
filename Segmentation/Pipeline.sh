#!/bin/bash 

cd $dex/programs/Segmentation; 

qcwd -l mem_free=22G,h_vmem=23G -N INDIV Run_Indiv_Model.sh

### Maxvmem = 31G
qcwd -l mem_free=30G,h_vmem=31G -N AGG -hold_jid INDIV Aggregate_Model.sh

qcwd -l mem_free=10G,h_vmem=12G -N COLL -hold_jid INDIV -hold_jid AGG Collapse_Predictor_Models.sh

qcwd -l mem_free=30G,h_vmem=31G -N CUT -hold_jid COLL Get_Model_Cutoffs.sh

qcwd -l mem_free=20G,h_vmem=21G -N PRED -hold_jid CUT Predict_From_Model.sh

Rbatch Get_Results.R -l mem_free=8G,h_vmem=10G -N RES -hold_jid PRED


