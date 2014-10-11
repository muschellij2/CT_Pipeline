rundir = pwd;
updir = fileparts(rundir);
outdir = fullfile(updir, 'Converted');
addon = '';
newsyntax = '1.2.840.10008.1.2';
copmode = 1;
createmode = 'Copy'; % could be 'Create';
x = dir(rundir);
nodir = ~[x.isdir]';
x = x(nodir);
x = {x.name}';
ifile = 1;
for ifile = 1:length(x)
    name = x{ifile};
    stub = regexprep(name, '(.*)[.].*', '$1');
    stub = [stub addon '.dcm']; %#ok<AGROW>
    name = fullfile(rundir, name);
    stub = fullfile(outdir, stub);
    dinfo = dicominfo(name);
    dinfo.TransferSyntaxUID = newsyntax;
    X = dicomread(name);
    dicomwrite(X, stub, dinfo, 'CreateMode', createmode);
end