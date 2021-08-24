function [UN] = getAvgUtilChain(self,U)
% [UN] = GETAVGUTILCHAIN(SELF,U)
% Return average utilization aggregated by chain
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

qn = self.model.getStruct();
%if nargin == 1
%    [Q] = self.model.getAvgHandles();
%end
[UNclass] = self.getAvgUtil();

% compute average chain metrics
UN = zeros(qn.nstations, qn.nchains);
for c=1:qn.nchains
    inchain = find(qn.chains(c,:));
    for i=1:qn.nstations
        if ~isempty(UNclass)
            UN(i,c) = sum(UNclass(i,inchain)); %#ok<FNDSB>
        end
    end
end
end
