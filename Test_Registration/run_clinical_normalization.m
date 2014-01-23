clear;
addpath('/home/bst/student/jmuschel')
pathdef;
cd('/dexter/disk2/smart/stroke_ct/ident/programs/Test_Registration');
spm_jobman('initcfg');

spm_jobman('run','CT_Normalize_All_BB.mat');
exit