function [qnchains,chains] = getChains(self, rt)
% [QNCHAINS,CHAINS] = GETCHAINS(RT)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

if ~exist('rt','var')
    rt = self.getRoutingMatrix();
end
refstat = self.getReferenceStations();
% getChains
[C, inChain] = weaklyconncomp(rt+rt');
chainCandidates = cell(1,C);
for c=1:C
    chainCandidates{c} = find(inChain==c);
end

K = getNumberOfClasses(self);
M = getNumberOfStations(self);

chains = [];
for t=1:length(chainCandidates)
    if length(chainCandidates{t})>1
        chains(end+1,unique(mod(chainCandidates{t}-1,K)+1)) = 1;
    end
end
try
    chains = sortrows(chains,'descend');
catch
    chains = sortrows(chains,-(1:size(chains,2)));
end

for c=1:size(chains,1)
    inchain = find(chains(c,:));
    if sum(refstat(inchain) == refstat(inchain(1))) ~= length(inchain)
        refstat(inchain) = refstat(inchain(1));
        %        line_error(sprintf('Classes in chain %d have different reference stations. Chain %d classes: %s', c, c, int2str(inchain)));
    end
end

visits = cell(size(chains,1),1); % visits{c}(i,j) is the number of visits that a chain-c job pays at node i in class j
for c=1:size(chains,1)
    inchain = find(chains(c,:));
    cols = [];
    for i=1:M
        for k=inchain(:)'
            cols(end+1) = (i-1)*K+k;
        end
    end
    Pchain = rt(cols,cols); % routing probability of the chain
    visited = nansum(Pchain,2) > 0;
    
    %                Pchain(visited,visited)
    %                if ~dtmc_isfeasible(Pchain(visited,visited))
    %                    line_error(sprintf('The routing matrix in chain %d is not stochastic. Chain %d classes: %s',c, c, int2str(inchain)));
    %                end
    alpha_visited = dtmc_solve(Pchain(visited,visited));
    alpha = zeros(1,M*K); alpha(visited) = alpha_visited;
    if max(alpha)>=1-1e-10
        line_error(mfilename,'One chain has an absorbing state.');
    end
    
    visits{c} = zeros(M,K);
    for i=1:M
        for k=1:length(inchain)
            visits{c}(i,inchain(k)) = alpha((i-1)*length(inchain)+k);
        end
    end
    visits{c} = visits{c} / sum(visits{c}(refstat(inchain(1)),inchain));
end
qnchains = cell(1, size(chains,1));
for c=1:size(chains,1)
    qnchains{c} = Chain(['Chain',num2str(c)]);
    for r=find(chains(c,:))
        qnchains{c}.addClass(self.classes{r}, visits{c}(:,r), r);
    end
end

end
