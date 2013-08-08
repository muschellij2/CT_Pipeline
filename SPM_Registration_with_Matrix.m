cd('~/CT_Registration/ICES/205-517/Separated/SPM/');
refname = '205-517_20110215_1701_Added.nii';
volname = '205-517_20110215_2230_2_Added.nii';
ref = spm_vol(refname);
vol = spm_vol(volname);

x = spm_coreg(ref, vol);
M = spm_matrix(x);
spm_get_space(volname, M * vol.mat);

% for re-slicing i usually construct my flags first:
flags= struct('interp',5,'mask',1,'mean',0,'which',1,'wrap',[0 0 0]');
% then the files i want to re-slice, source file first (in our example T1 image)
files = {refname; volname};
spm_reslice(files, flags);
