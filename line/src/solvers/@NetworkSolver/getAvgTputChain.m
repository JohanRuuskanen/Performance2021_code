function [TN] = getAvgTputChain(self,T)
% [TN] = GETAVGTPUTCHAIN(SELF,T)

% Return average throughputs aggregated by chain
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

qn = self.model.getStruct();
%if nargin == 1
%    [Q] = self.model.getAvgHandles();
%end
[TNclass] = self.getAvgTput();

% compute average chain metrics
TN = zeros(qn.nstations, qn.nchains);
for c=1:qn.nchains
    inchain = find(qn.chains(c,:));
    for i=1:qn.nstations
        if ~isempty(TNclass)
            TN(i,c) = sum(TNclass(i,inchain)); %#ok<FNDSB>
        end
    end
end
end
