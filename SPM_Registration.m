clear;
% List of open inputs
% Coreg: Estimate & Reslice: Reference Image - cfg_files
% Coreg: Estimate & Reslice: Source Image - cfg_files
basedir='/Users/muschellij2/Dropbox/CTR/DHanley/CT_Registration/ICES';
jobfile = {fullfile(basedir, 'programs', 'SPM_Registration_job.m')};
load(fullfile(basedir, 'Scans.mat'));
nrun = size(res.patientName, 1); % enter the number of runs here

% jobs = repmat(jobfile, 1, nrun);
% inputs = cell(2, nrun);

spm('defaults', 'FMRI');

for crun = 5:nrun
    rundir=fileparts(res.PostOp{crun});
    cd(rundir);
    jobs = repmat(jobfile, 1, 1);
    inputs = cell(2, 1);
    inputs{1, 1} = res.PreOp(crun); % Coreg: Estimate & Reslice: Reference Image - cfg_files
    inputs{2, 1} = res.PostOp(crun); % Coreg: Estimate & Reslice: Source Image - cfg_files
%     inputs{1, crun} = res.PreOp(crun); % Coreg: Estimate & Reslice: Reference Image - cfg_files
%     inputs{2, crun} = res.PostOp(crun); % Coreg: Estimate & Reslice: Source Image - cfg_files
    spm_jobman('serial', jobs, '', inputs{:});
    cmd=sprintf('gzip "%s"/*.nii',  rundir);
    system(cmd);
    cmd=sprintf('rm "%s"/*.nii',  rundir);
    system(cmd);    
end
