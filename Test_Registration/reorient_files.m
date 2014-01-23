

%[P, sts] = spm_select([1 Inf], 'image','Images to reorient');
%if ~sts, return; else P = cellstr(P); end

redir = '/dexter/disk2/smart/stroke_ct/ident/Test_Registration/reoriented/';
cd(redir);
f = dir(redir);
f = regexpi({f.name},'^\d\d.*\.nii|','match');
P = [f{:}];
for e = 1:numel(P)
    x = P{e};
    P{e} = fullfile(redir, x);
end

for i=1:numel(P)
    Mats(:,:,i) = spm_get_space(P{i});
end
for i=1:numel(P)
    fname=[P{i}(1:(end-4)),'_reorient.mat'];
    x = load(fname, '-ascii');
    mat = spm_matrix(x);
    spm_get_space(P{i},mat*Mats(:,:,i));
    disp(P{i});
end