
% clear;
% addpath('/home/bst/student/jmuschel')
% pathdef;
dexpath='/dexter/disk2/smart/stroke_ct/ident/programs/Test_Registration'
spm_jobman('initcfg');
matfile = fullfile(dexpath, 'CT_Normalize_All_BB.mat')
spm_jobman('run', matfile);
exit