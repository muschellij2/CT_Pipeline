
% clear;
% addpath('/home/bst/student/jmuschel')
% pathdef;
clear;
dexpath='/dexter/disk2/smart/stroke_ct/ident/programs/Test_Registration';
outdir = fullfile(dexpath, 'reoriented');

matfile = fullfile(dexpath, 'CT_Normalize_All_BB.mat')
l = load(matfile);
imgs = [l.matlabbatch{1}.spm.tools.MRI.CTnorm.images];
lesions = [l.matlabbatch{1}.spm.tools.MRI.CTnorm.ctles];
iimg  =1;
for iimg = 1:length(imgs)
	runimgs = strvcat(imgs{iimg}, lesions{iimg});
	nii_setorigin(imgs{iimg});
	% disp(imgs{iimg});
end

spm_jobman('initcfg');
spm_jobman('run', matfile);
exit