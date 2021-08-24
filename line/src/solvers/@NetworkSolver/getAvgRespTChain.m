function [RN] = getAvgRespTChain(self,R)
% [RN] = GETAVGRESPTCHAIN(SELF,R)
% Return average response time aggregated by chain
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

qn = self.model.getStruct();
%if nargin == 1
%    [Q] = self.model.getAvgHandles();
%end
[RNclass] = self.getAvgRespT();

% compute chain visits
alpha = zeros(qn.nstations,qn.nclasses);
for c=1:qn.nchains
    inchain = find(qn.chains(c,:));
    for i=1:qn.nstations
        for k=inchain % for all classes within the chain ( a class belongs to a single chain, the reference station must be identical
            %                        for all classes within a chain )
            alpha(i,k) = alpha(i,k) + qn.visits{c}(i,k)/sum(qn.visits{c}(qn.refstat(k),inchain));
        end
    end
end
alpha(~isfinite(alpha))=0;

% compute average chain metrics
RN = zeros(qn.nstations, qn.nchains);
for c=1:qn.nchains
    inchain = find(qn.chains(c,:));
    for i=1:qn.nstations
        if ~isempty(RNclass)
            RN(i,c) = RNclass(i,inchain)*alpha(i,inchain)';
        end
    end
end
end
