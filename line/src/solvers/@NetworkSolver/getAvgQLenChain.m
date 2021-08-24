function [QN] = getAvgQLenChain(self,Q)
% [QN] = GETAVGQLENCHAIN(SELF,Q)

% Return average queue-lengths aggregated by chain
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

qn = self.model.getStruct();
%if nargin == 1
%    [Q] = self.model.getAvgHandles();
%end
[QNclass] = self.getAvgQLen();

% compute average chain metrics
QN = zeros(qn.nstations, qn.nchains);
for c=1:qn.nchains
    inchain = find(qn.chains(c,:));
    for i=1:qn.nstations
        if ~isempty(QNclass)
            QN(i,c) = sum(QNclass(i,inchain));
        end
    end
end
end
