function model = PMIF2LINE(filename,modelName)
% MODEL = PMIF2LINE(FILENAME,MODELNAME)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
verbose = false;
qn = PMIF2QN(filename,verbose);
model = QN2LINE(qn, modelName);
end
