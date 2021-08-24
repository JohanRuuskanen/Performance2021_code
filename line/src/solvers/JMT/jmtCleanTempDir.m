pd = pwd;
cwd = fullfile(tempdir,'jsimg');
cd(cwd)
delete *.jsimg
delete *.jsimg-result.jsim
cwd = fullfile(tempdir,'jmva');
cd(cwd)
delete *.jmva
delete *.jsimg-result.jsim
cd(pd)
