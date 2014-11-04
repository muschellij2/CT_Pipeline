#!/bin/bash 

cd $dex/programs/Segmentation; 

Rbatch Create_Filename_Rda.R -N FILES -l mem_free=2G,h_vmem=3G

Rbatch Sample_and_Result_Format.R -l mem_free=2G,h_vmem=3G \
-N FORMATS -hold_jid FILES

#############################
# Get true volume from original data predictor matrix
#############################
Rbatch 'Create_Predictors.R -c none -r FALSE -o TRUE' \
-hold_jid FILES \
-l mem_free=30G,h_vmem=33G -N CREATE_ORIG -t 1-112


Rbatch Get_True_Vol.R -N TRUTH -hold_jid CREATE_ORIG \
-l mem_free=16G,h_vmem=17G

###############################
# Run N3 Correction
###############################
# qcwd -l mem_free=12G,h_vmem=13G -hold_jid FILES \
# -N N3 N3_Correct_Images.sh

Rbatch 'Create_Predictors.R -c N3_SS -r FALSE -o TRUE' \
-hold_jid N3 -l mem_free=30G,h_vmem=33G -N CREATE_N3 -t 1-112

Rbatch 'Create_Predictors.R -c N4_SS -r FALSE -o TRUE' \
-hold_jid FILES -l mem_free=30G,h_vmem=33G -N CREATE_N4 -t 1-112

###############################
# Run Registrations
###############################
qcwd -l mem_free=12G,h_vmem=13G -N REG -hold_jid FILES \
CT_Register_to_Template_SyN.sh

###############################
# Run Thresholding determination and threshold
###############################
Rbatch Determine_Threshold_Registration_Values.R -hold_jid REG \
-N GET_REG_THRESH -l mem_free=3G,h_vmem=4G

Rbatch Threshold_Registration_Values.R -l mem_free=8G,h_vmem=10G \
-N THRESH -hold_jid GET_REG_THRESH -t 1-112

# #######################################
# # Need to figure out Registration threshold here
# #######################################
Rbatch 'Create_Predictors.R -c Rigid -r TRUE -o TRUE' \
-hold_jid_ad THRESH \
-l mem_free=30G,h_vmem=33G -N CREATE_Rigid -t 1-112

Rbatch 'Create_Predictors.R -c Rigid_sinc -r TRUE -o TRUE' \
-hold_jid_ad THRESH \
-l mem_free=30G,h_vmem=33G -N CREATE_Rigid_sinc -t 1-112


qcwd -l mem_free=30G,h_vmem=33G -N INDIV \
-hold_jid_ad CREATE_N3 -hold_jid_ad CREATE_N4 \
-hold_jid_ad CREATE_Rigid -hold_jid_ad CREATE_Rigid_sinc \
-hold_jid_ad CREATE_ORIG Run_Indiv_Model.sh



# ### Maxvmem = 31G
qcwd -l mem_free=30G,h_vmem=31G -N AGGDATA \
-hold_jid FORMATS -hold_jid CREATE_Rigid_sinc \
Create_Aggregate_Data.sh


qcwd -l mem_free=10G,h_vmem=12G -N INCLUDE \
-hold_jid INDIV -hold_jid_ad AGGDATA Inclusion_Criteria_Results.sh

qcwd -l mem_free=15G,h_vmem=17G -N PLOT \
-hold_jid AGGDATA Plot_Agg_Relationship.sh


qcwd -l mem_free=30G,h_vmem=31G -N AGG -hold_jid INDIV \
-hold_jid AGGDATA Aggregate_Model.sh

qcwd -l mem_free=10G,h_vmem=12G -N COLL -hold_jid FORMATS \
-hold_jid AGG Collapse_Predictor_Models.sh



qcwd -l mem_free=30G,h_vmem=31G -N CUT \
-hold_jid COLL Get_Model_Cutoffs.sh

qcwd -l mem_free=24G,h_vmem=25G -N SCUT \
-hold_jid COLL -hold_jid_ad CUT Get_Smooth_Model_Cutoffs.sh

qcwd -l mem_free=20G,h_vmem=21G -N PRED \
-hold_jid CUT -hold_jid SCUT Predict_From_Model.sh

# Rbatch Get_Results.R -l mem_free=8G,h_vmem=10G \
# -N RES -hold_jid PRED


