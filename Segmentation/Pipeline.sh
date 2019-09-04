#!/bin/bash 

cd $dex/programs/Segmentation; 

Rbatch Create_Filename_Rda.R -N FILES -l mem_free=2G,h_vmem=3G

Rbatch Sample_and_Result_Format.R -l mem_free=2G,h_vmem=3G \
-N FORMATS -hold_jid FILES

qcwd -l mem_free=20G,h_vmem=21G -N ZSCORE \
-hold_jid FORMATS \
Zscore_Template.sh

################################
# Create Flipped Image
################################
Rbatch 'Flipped_Stroke_Images.R -c none -r TRUE' \
-hold_jid_ad FORMATS \
-l mem_free=10G,h_vmem=13G -N FLIP -t 1-112

################################
Rbatch 'Flipped_Stroke_Images.R -c Rigid -r TRUE' \
-hold_jid_ad FORMATS \
-l mem_free=10G,h_vmem=13G -N FLIP_Rigid -t 1-112


#############################
# Get true volume from original data predictor matrix
#############################
Rbatch 'Create_Predictors.R -c none -r TRUE -o TRUE' \
-hold_jid_ad ZSCORE \
-l mem_free=30G,h_vmem=33G -N CREATE_ORIG -t 1-112
# qdel 4343869.111

Rbatch 'Create_Predictors.R -c none -r TRUE -o TRUE' \
-hold_jid ZSCORE \
-l mem_free=60G,h_vmem=62G -N CREATE_ORIG_111 -t 111


Rbatch Get_True_Vol.R -N TRUTH -hold_jid CREATE_ORIG \
-l mem_free=16G,h_vmem=17G

###############################
# Run N3 Correction
###############################
# qcwd -l mem_free=12G,h_vmem=13G -hold_jid FILES \
# -N N3 N3_Correct_Images.sh

# Rbatch 'Create_Predictors.R -c N3_SS -r FALSE -o TRUE' \
# -hold_jid N3 -l mem_free=30G,h_vmem=33G \
# -N CREATE_N3 -t 1-112

# Rbatch 'Create_Predictors.R -c N4_SS -r TRUE -o TRUE' \
# -hold_jid FILES -l mem_free=30G,h_vmem=33G -N \
# CREATE_N4 -t 1-112

###############################
# Run Registrations
###############################
qcwd -l mem_free=12G,h_vmem=13G -N REG -hold_jid_ad ZSCORE \
CT_Register_to_Template_SyN.sh

# ###############################
# # Run Thresholding determination and threshold
# ###############################
# Rbatch Determine_Threshold_Registration_Values.R \
# -hold_jid REG \
# -N GET_REG_THRESH -l mem_free=5G,h_vmem=6G

# Rbatch Threshold_Registration_Values.R \
# -l mem_free=8G,h_vmem=10G \
# -N THRESH -hold_jid GET_REG_THRESH -t 1-112

# #######################################
# # Need to figure out Registration threshold here
# #######################################
Rbatch 'Create_Predictors.R -c Rigid -r TRUE -o TRUE' \
-hold_jid_ad THRESH -hold_jid_ad REG \
-l mem_free=30G,h_vmem=33G -N CREATE_Rigid -t 1-112

# Rbatch 'Create_Predictors.R -c Rigid_sinc -r FALSE -o TRUE' \
# -hold_jid_ad THRESH -hold_jid_ad REG \
# -l mem_free=30G,h_vmem=33G -N CREATE_Rigid_sinc -t 1-112


#############################
# Rerun winsorized if you forget again
#############################
# Rbatch 'Rerun_Winsorized.R -c none -r TRUE -o TRUE' \
# -hold_jid_ad CREATE_ORIG \
# -l mem_free=10G,h_vmem=12G -N RERUN -t 1-112

# Rbatch 'Rerun_Winsorized.R -c Rigid -r TRUE -o TRUE' \
# -hold_jid_ad CREATE_RIGID \
# -l mem_free=5G,h_vmem=8G -N RERUN_RIGID -t 1-112

# Rbatch 'Rerun_Winsorized.R -c Rigid_sinc -r TRUE -o TRUE' \
# -hold_jid_ad CREATE_RIGID_sinc \
# -l mem_free=5G,h_vmem=8G -N RERUN_SRIGID -t 1-112

###############################
# Creating Aggregate Data frame from training
###############################
# ### Maxvmem = 31G
qcwd -l mem_free=30G,h_vmem=31G -N AGGDATA \
-hold_jid FORMATS -hold_jid CREATE_Rigid \
-hold_jid CREATE_ORIG \
Create_Aggregate_Data.sh

###############################
# Creating Inclusion/Exclusion from distribution of data
###############################
qcwd -l mem_free=10G,h_vmem=12G -N INCLUDE \
 -hold_jid_ad AGGDATA Inclusion_Criteria_Results.sh
 
###############################
# Creating Validation Data frame from training
###############################
qcwd -l mem_free=30G,h_vmem=31G -N VALIDDATA \
-hold_jid FORMATS -hold_jid CREATE_Rigid_sinc \
-hold_jid CREATE_Rigid -hold_jid CREATE_ORIG \
-hold_jid_ad AGGDATA Create_Validation_Data.sh


###############################
# Creating Individual Models
###############################
# qcwd -l mem_free=30G,h_vmem=33G -N INDIV \
# -hold_jid_ad CREATE_N3 -hold_jid_ad CREATE_N4 \
# -hold_jid_ad CREATE_Rigid -hold_jid_ad CREATE_Rigid_sinc \
# -hold_jid_ad CREATE_ORIG -hold_jid AGGDATA \
# Run_Indiv_Model.sh

Rbatch 'Run_Indiv_Model.R -c Rigid' \
-l mem_free=10G,h_vmem=13G -N INDIV_Rig \
-hold_jid_ad CREATE_Rigid -hold_jid AGGDATA -t 1-112

# Rbatch 'Run_Indiv_Model.R -c Rigid_sinc' \
# -l mem_free=10G,h_vmem=13G -N INDIV_Rig_sinc \
# -hold_jid_ad CREATE_Rigid_sinc -hold_jid AGGDATA  -t 1-112

Rbatch 'Run_Indiv_Model.R -c none' \
-l mem_free=30G,h_vmem=33G -N INDIV \
-hold_jid_ad CREATE_ORIG -hold_jid AGGDATA -t 1-112

# Rbatch 'Run_Indiv_Model.R -c none' \
# -l mem_free=30G,h_vmem=33G -N VINDIV \
# -hold_jid_ad CREATE_ORIG -hold_jid AGGDATA  -t 11-112

# Rbatch 'Run_Indiv_Model.R -c N3_SS' \
# -l mem_free=30G,h_vmem=33G -N INDIV_N3 \
# -hold_jid_ad CREATE_N3 -hold_jid AGGDATA  -t 1-112

# Rbatch 'Run_Indiv_Model.R -c N4_SS' \
# -l mem_free=30G,h_vmem=33G -N INDIV_N4 \
# -hold_jid_ad CREATE_N4 -hold_jid AGGDATA  -t 1-112


Rbatch 'Get_SS_Statistics.R' \
-N STATS -l mem_free=10G,h_vmem=14G \
-hold_jid CUT -hold_jid SCUT -t 111

# ******STOPPED HERE*******

###############################
# Visualize relationship of data
###############################
qcwd -l mem_free=15G,h_vmem=17G -N PLOT \
-hold_jid AGGDATA Plot_Agg_Relationship.sh


###############################
# Create Aggregrate Model
###############################
qcwd -l mem_free=60G,h_vmem=61G -N AGG \
-hold_jid AGGDATA Aggregate_Model.sh


qcwd -l mem_free=10G,h_vmem=12G -N COLL \
-hold_jid FORMATS \
-hold_jid_ad AGG Collapse_Predictor_Models.sh

qcwd -l mem_free=10G,h_vmem=12G -N CUT \
-hold_jid_ad COLL -hold_jid INDIV Get_Model_Cutoffs.sh

qcwd -l mem_free=24G,h_vmem=25G -N SCUT \
-hold_jid_ad COLL -hold_jid_ad CUT \
Get_Smooth_Model_Cutoffs.sh


Rbatch 'Predict_From_Model.R -c Rigid' \
-N PRED_Rig -l mem_free=23G,h_vmem=25G \
-hold_jid CUT -hold_jid SCUT -t 1-112

# Rbatch 'Predict_From_Model.R -c Rigid_sinc' \
# -N PRED_Rig_sinc -l mem_free=20G,h_vmem=31G \
# -hold_jid CUT -hold_jid SCUT -t 1-112

Rbatch 'Predict_From_Model.R -c none' \
-N PRED -l mem_free=30G,h_vmem=32G \
-hold_jid CUT -hold_jid SCUT -t 1-112

Rbatch 'Predict_From_Model.R -c none' \
-N PRED_111 -l mem_free=69G,h_vmem=70G \
-hold_jid CUT -hold_jid SCUT -t 111



# Rbatch 'Predict_From_Model.R -c N3_SS' \
# -N PRED_N3 -l mem_free=20G,h_vmem=31G \
# -hold_jid CUT -hold_jid SCUT -t 1-112

# Rbatch 'Predict_From_Model.R -c N4_SS' \
# -N PRED_N4 -l mem_free=20G,h_vmem=31G \
# -hold_jid CUT -hold_jid SCUT -t 1-112

# qcwd -l mem_free=20G,h_vmem=31G -N PRED \
# -hold_jid CUT -hold_jid SCUT Predict_From_Model.sh

# i=Rigid;
for i in Rigid Rigid_sinc none N3_SS N4; do 
	x="Select_Aggregate_Model.R -c ${i}";
	Rbatch "${x}" \
	-N SELECT_${i} -l mem_free=4G,h_vmem=5G;
done;

for i in Rigid Rigid_sinc none N3_SS N4; do 
	x="Select_Aggregate_Model_Results.R -c ${i}";
	Rbatch "${x}" \
	-N SRESULT_${i} -l mem_free=30G,h_vmem=31G \
	-hold_jid SELECT_${i} -t 1-55;	
done;

Rbatch Get_Results.R -l mem_free=8G,h_vmem=10G \
-N RES -hold_jid PRED -hold_jid PRED_N3 -hold_jid PRED_N4 \
-hold_jid PRED_Rig -hold_jid PRED_Rig_sinc


####### need final predictor
Rbatch Create_Final_Predictor_Image.R -N FINAL \
-l mem_free=7G,h_vmem=8G \
-hold_jid PRED_Rig -t 1-6;

Rbatch Create_Final_Predictor_Image.R -N FINAL \
-l mem_free=7G,h_vmem=8G \
-hold_jid PRED -t 7-12;

qcwd -l mem_free=5G,h_vmem=6G -N TRANS \
-hold_jid FILES \
-hold_jid FINAL Transform_Masks_Back.sh

# Rbatch Final_Predictor_Volumes_Native.R -N FVOLS \
# -l mem_free=30G,h_vmem=33G -hold_jid TRANS

Rnosave Final_Predictor_Cutoff_Volumes_Native.R -N FCUTVOLS \
-l mem_free=18G,h_vmem=20G -hold_jid TRANS

Rnosave Native_Prediction_Volumes.R -N NATIVEVOLS \
-l mem_free=8G,h_vmem=12G -hold_jid TRANS

qcwd -N RESFINAL -l mem_free=20G,h_vmem=23G -hold_jid TRANS \
Create_Results_Final.sh 


Rnosave Figures_for_Paper.R -N FIGURES \
-l mem_free=8G,h_vmem=12G -hold_jid RESFINAL
