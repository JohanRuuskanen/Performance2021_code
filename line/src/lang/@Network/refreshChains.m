function [chains, visits, rt] = refreshChains(self, wantVisits, wantCapacity)
% [CHAINS, VISITS, RT] = REFRESHCHAINS(WANTVISITS, WANTCAPACITY)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

if ~exist('wantVisits','var')
    wantVisits = true;
end

if nargin == 1
    if isempty(self.qn)
        line_error(mfilename,'refreshRoutingMatrix cannot retrieve station rates, please pass them as an input parameters.');
    end
end

rates = self.qn.rates;
I = self.getNumberOfNodes();
M = self.getNumberOfStatefulNodes();
K = self.getNumberOfClasses();
refstat = self.getReferenceStations();
[rt,~,csmask, rtnodes] = self.refreshRoutingMatrix(rates);
% getChains
[C,inChain] = weaklyconncomp(csmask+csmask');
%[C,inChain] = weaklyconncomp(rt+rt')

chainCandidates = cell(1,C);
for c=1:C
    chainCandidates{c} = find(inChain==c);
end

chains = [];
for t=1:length(chainCandidates)
    %    if length(chainCandidates{t})>1
    %        chains(end+1,unique(mod(chainCandidates{t}-1,K)+1)) = true;
    chains(end+1,chainCandidates{t}) = true;
    %    end
end

try
    chains = sortrows(chains,'descend');
catch
    chains = sortrows(chains);
end

for c=1:size(chains,1)
    inchain{c} = find(chains(c,:));
end

for c=1:size(chains,1)
    if sum(refstat(inchain{c}) == refstat(inchain{c}(1))) ~= length(inchain{c})
        refstat(inchain{c}) = refstat(inchain{c}(1));
        %        line_error(sprintf('Classes in chain %d have different reference stations. Chain %d classes: %s', c, c, int2str(inchain{c})));
    end
end

%% generate station visits
visits = cell(size(chains,1),1); % visits{c}(i,j) is the number of visits that a chain-c job pays at node i in class j
if wantVisits
    for c=1:size(chains,1)
        cols = [];
        for i=1:M
            for k=inchain{c}(:)'
                cols(end+1) = (i-1)*K+k;
            end
        end
        Pchain = rt(cols,cols); % routing probability of the chain
        visited = sum(Pchain,2) > 0;
        %                Pchain(visited,visited)
        %                if ~dtmc_isfeasible(Pchain(visited,visited))
        %                    line_error(sprintf('The routing matrix in chain %d is not stochastic. Chain %d classes: %s',c, c, int2str(inchain{c})));
        %                end
        alpha_visited = dtmc_solve(Pchain(visited,visited));
        alpha = zeros(1,M*K); alpha(visited) = alpha_visited;
        if max(alpha)>=1-1e-10
            %disabled because a self-looping customer is an absorbing chain
            %line_error(mfilename,'Line:ChainAbsorbingState','One chain has an absorbing state.');
        end
        visits{c} = zeros(M,K);
        for i=1:M
            for k=1:length(inchain{c})
                visits{c}(i,inchain{c}(k)) = alpha((i-1)*length(inchain{c})+k);
            end
        end
        visits{c} = visits{c} / sum(visits{c}(refstat(inchain{c}(1)),inchain{c}));
        visits{c} = abs(visits{c});
    end
end

%% generate node visits
nchains = size(chains,1);
if wantVisits
    nodevisits = cell(1,nchains);
    for c=1:nchains
        nodes_cols = [];
        for i=1:I
            for k=inchain{c}(:)'
                nodes_cols(end+1) = (i-1)*K+k;
            end
        end
        nodes_Pchain = rtnodes(nodes_cols, nodes_cols); % routing probability of the chain
        nodes_visited = sum(nodes_Pchain,2) > 0;
        
        nodes_alpha_visited = dtmc_solve(nodes_Pchain(nodes_visited,nodes_visited));
        nodes_alpha = zeros(1,I); nodes_alpha(nodes_visited) = nodes_alpha_visited;
        nodevisits{c} = zeros(I,K);
        for i=1:I
            for k=1:length(inchain{c})
                nodevisits{c}(i,inchain{c}(k)) = nodes_alpha((i-1)*length(inchain{c})+k);
            end
        end
        nodevisits{c} = nodevisits{c} / sum(nodevisits{c}(refstat(inchain{c}(1)),inchain{c}));
        nodevisits{c}(nodevisits{c}<0) = 0; % remove small numerical perturbations
    end
    self.qn.nodevisits = nodevisits;
end

if ~isempty(self.qn) %&& isprop(self.qn,'chains')
    if exist('nodevisits','var')
        self.qn.setChains(chains, visits, rt, nodevisits);
    else
        self.qn.setChains(chains, visits, rt , {});
    end
end

%% call dependent refreshes
self.refreshCapacity(); % capacity depends on chains and rates

for c=1:self.qn.nchains
    self.qn.visits{c}(isnan(self.qn.visits{c})) = 0;
end

for r=1:self.qn.nclasses
    if all(self.qn.routing(:,r) == -1)
        line_error(sprintf('Routing strategy in class %d is unspecified at all nodes.',r));
    end
end

end