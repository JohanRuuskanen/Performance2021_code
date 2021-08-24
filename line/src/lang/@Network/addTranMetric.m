function addTranMetric(self, perfIndex)
% ADDTRANMETRIC(PERFINDEX)
% Add handle for a transient metric

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
nel = numel(perfIndex);
if nel == 1
    self.perfIndex.Tran{end+1,1} = perfIndex;
else
    self.perfIndex.Tran(end+1:end+numel(perfIndex), 1) = perfIndex(:);
end
end
