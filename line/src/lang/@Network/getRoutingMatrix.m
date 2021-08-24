function [rt,rtNodes,connMatrix,rtNodesByClass,rtNodesByStation] = getRoutingMatrix(self, arvRates)
% [RT,RTNODES,CONNMATRIX,RTNODESBYCLASS,RTNODESBYSTATION] = GETROUTINGMATRIX(ARVRATES)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

if ~exist('arvRates','var')
    for r=self.getIndexOpenClasses
        arvRates(r) = 1 / self.getSource.input.sourceClasses{r}{end}.getMean;
    end
end

nodeNames = self.getNodeNames();
% connectivity matrix
connMatrix = zeros(self.getNumberOfNodes);
for r=1:size(self.links,1)
    connMatrix(self.links{r}{1}, self.links{r}{2}) = 1;
end
rtNodesByClass = {};
rtNodesByStation = {};
hasOpenClasses = self.hasOpenClasses();
if ~exist('arvRates','var')
    if hasOpenClasses
        line_error(mfilename,'getRoutingMatrix requires arrival rates for open classes.');
    end
end

M = self.getNumberOfNodes;
K = self.getNumberOfClasses;
NK = self.getNumberOfJobs;
rtNodes = zeros(M*K);
% The first loop considers the class at which a job enters the
% target station
for i=1:M
    switch class(self.nodes{i}.output)
        case 'Forker'
            for j=1:M
                for k=1:K
                    if connMatrix(i,j)>0
                        rtNodes((i-1)*K+k,(j-1)*K+k)=1.0;
                        switch self.nodes{i}.output.outputStrategy{k}{2}
                            case RoutingStrategy.PROB
                                if length(self.nodes{i}.output.outputStrategy{k}{end}) ~= sum(connMatrix(i,:))
                                    line_error(mfilename,'Fork must have 1.0 routing probability towards all outgoing links.');
                                end
                                for t=1:length(self.nodes{i}.output.outputStrategy{k}{end}) % for all outgoing links
                                    if self.nodes{i}.output.outputStrategy{k}{end}{t}{2} ~= 1.0
                                        line_error(mfilename,'Fork must have 1.0 routing probability towards all outgoing links, but a routing probability is at %f.',self.nodes{i}.output.outputStrategy{k}{end}{t}{2});
                                    end
                                end
                        end
                    end
                end
            end
        otherwise
            for k=1:K
                switch self.nodes{i}.output.outputStrategy{k}{2}
                    case RoutingStrategy.PROB
                        if isinf(NK(k)) || ~isa(self.nodes{i},'Sink') 
                            %rtNodes((i-1)*K+k,(j-1)*K+k) = self.modifiedRoutingTable{k,k}(i,j);
                            for t=1:length(self.nodes{i}.output.outputStrategy{k}{end}) % for all outgoing links
                                j = findstring(nodeNames, self.nodes{i}.output.outputStrategy{k}{end}{t}{1}.name);
                                rtNodes((i-1)*K+k,(j-1)*K+k) = self.nodes{i}.output.outputStrategy{k}{end}{t}{2};
                            end
                        end
                    case RoutingStrategy.DISABLED
                        % set a small numerical tolerance to avoid messing
                        % up with routing chain CTMC solution
                        for j=1:M
                            if connMatrix(i,j)>0
                                %rtNodes((i-1)*K+k,(j-1)*K+k) = Distrib.Zero;
                            end
                        end
                    case {RoutingStrategy.RAND, RoutingStrategy.RRB, RoutingStrategy.JSQ}
                        if isinf(NK(k)) % open class
                            for j=1:M
                                if connMatrix(i,j)>0
                                    rtNodes((i-1)*K+k,(j-1)*K+k)=1/sum(connMatrix(i,:));
                                end
                            end
                        elseif ~isa(self.nodes{i},'Source') && ~isa(self.nodes{i},'Sink') % don't route closed classes out of source nodes
                            connMatrixClosed = connMatrix;
                            if connMatrixClosed(i,self.getNodeIndex(self.getSink))
                                connMatrixClosed(i,self.getNodeIndex(self.getSink)) = 0;
                            end
                            for j=1:M
                                if connMatrixClosed(i,j)>0
                                    rtNodes((i-1)*K+k,(j-1)*K+k)=1/(sum(connMatrixClosed(i,:)));
                                end
                            end
                        end
                    otherwise
                        for j=1:M
                            if connMatrix(i,j)>0
                                rtNodes((i-1)*K+k,(j-1)*K+k) = Distrib.Zero;
                            end
                        end
                        %line_error([self.nodes{i}.output.outputStrategy{k}{2},' routing policy is not yet supported.']);
                end
            end
    end
end

% The second loop corrects the first one at nodes that change
% the class of the job in the service section.


for i=1:self.getNumberOfNodes % source
    if isa(self.nodes{i}.server,'StatelessClassSwitcher')
        Pi = rtNodes(((i-1)*K+1):i*K,:);
        Pcs = self.nodes{i}.server.csFun(1:K,1:K);
        %for r=1:K
        %    for s=1:K
        %        Pcs(r,s) = self.nodes{i}.server.csFun(r,s);
        %    end
        %end
        rtNodes(((i-1)*K+1):i*K,:) = 0;
        for j=1:M % destination
            Pij = Pi(1:K,((j-1)*K+1):j*K); %Pij(r,s)            
           % for r=1:K
                for s=1:K
                    % Find the routing probability section determined by the router section in the first loop
                    %Pnodes(((i-1)*K+1):i*K,((j-1)*K+1):j*K) = Pcs*Pij;
                    rtNodes(((i-1)*K+1) : ((i-1)*K+K),(j-1)*K+s) = Pcs(1:K,s)*Pij(s,s);
                end
            %end
        end
    elseif isa(self.nodes{i}.server,'StatefulClassSwitcher')
        Pi = rtNodes(((i-1)*K+1):i*K,:);
        for r=1:K
            for s=1:K
                Pcs(r,s) = self.nodes{i}.server.csFun(r,s,[],[]); % get csmask
            end
        end
        rtNodes(((i-1)*K+1):i*K,:) = 0;
        if isa(self.nodes{i}.server,'CacheClassSwitcher')
            for r=1:K
                if (isempty(find(r == self.nodes{i}.server.hitClass)) && isempty(find(r == self.nodes{i}.server.missClass)))
                    Pcs(r,:) = Pcs(r,:)/sum(Pcs(r,:));
                end
            end
            
            for r=1:K
                if (isempty(find(r == self.nodes{i}.server.hitClass)) && isempty(find(r == self.nodes{i}.server.missClass)))
                    for j=1:M % destination
                        for s=1:K
                            Pi((i-1)*K+r,(j-1)*K+s) = 0;
                        end
                    end
                end
            end      
            
            for r=1:K
                if length(self.nodes{i}.server.actualHitProb)>=r && length(self.nodes{i}.server.hitClass)>=r
                    ph = self.nodes{i}.server.actualHitProb(r);
                    pm = self.nodes{i}.server.actualMissProb(r);
                    h = self.nodes{i}.server.hitClass(r);
                    m = self.nodes{i}.server.missClass(r);
                    rtNodes((i-1)*K+r,(i-1)*K+h) = ph;
                    rtNodes((i-1)*K+r,(i-1)*K+m) = pm;
                else
                    if length(self.nodes{i}.server.hitClass)>=r
                        h = self.nodes{i}.server.hitClass(r);
                        m = self.nodes{i}.server.missClass(r);
                        rtNodes((i-1)*K+r,(i-1)*K+h) = NaN;
                        rtNodes((i-1)*K+r,(i-1)*K+m) = NaN;
                    end
                end
            end
            
            for j=1:M % destination
                Pij = Pi(1:K,((j-1)*K+1):j*K); %Pij(r,s)
                for r=1:K
                    if ~(isempty(find(r == self.nodes{i}.server.hitClass)) && isempty(find(r == self.nodes{i}.server.missClass)))
                        for s=1:K
                            % Find the routing probability section determined by the router section in the first loop
                            %Pnodes(((i-1)*K+1):i*K,((j-1)*K+1):j*K) = Pcs*Pij;
                            rtNodes((i-1)*K+r,(j-1)*K+s) = Pcs(r,s)*Pij(s,s);
                        end
                    end
                end
            end
        end
    end
    
    % ignore all chains containing a Pnodes column that sums to 0,
    % since these are classes that cannot arrive to the node
    % unless this column belongs to the source
    colsToIgnore = find(sum(rtNodes,1)==0);
    if hasOpenClasses
        idxSource = self.getIndexSourceNode;
        colsToIgnore = setdiff(colsToIgnore,(idxSource-1)*K+(1:K));
    end
    
    % We route back from the sink to the source. Since open classes
    % have an infinite population, if there is a class switch QN
    % with the following chains
    % Source -> (A or B) -> C -> Sink
    % Source -> D -> Sink
    % We can re-route class C into the source either as A or B or C.
    % We here re-route back as C and leave for the chain analyzer
    % to detect that C is in a chain with A and B and change this
    % part.
        
    [C,inChain]=weaklyconncomp(rtNodes+rtNodes');
    inChain(colsToIgnore) = 0;
    chainCandidates = cell(1,C);
    for r=1:C
        chainCandidates{r} = find(inChain==r);
    end
    
    chainsPnodes = []; % columns are classes? rows are definitely chains
    for t=1:length(chainCandidates)
        if length(chainCandidates{t})>1
            chainsPnodes(end+1,unique(mod(chainCandidates{t}-1,K)+1))=1;
        end
    end
    try
        chainsPnodes = sortrows(chainsPnodes,'descend');
    catch % old MATLABs
        chainsPnodes = sortrows(chainsPnodes);
    end
    % this routes open classes back from the sink into the source
    % it will not work with non-renewal arrivals as it choses in which open
    % class to reroute a job with probability depending on the arrival rates
    if hasOpenClasses
        arvRates(isnan(arvRates)) = 0;
        idxSink = self.getIndexSinkNode;
        for s=self.getIndexOpenClasses
            s_chain = find(chainsPnodes(:,s));
            others_in_chain = find(chainsPnodes(s_chain,:));
            rtNodes((idxSink-1)*K+others_in_chain,(idxSource-1)*K+others_in_chain) = repmat(arvRates(others_in_chain)/sum(arvRates(others_in_chain)),length(others_in_chain),1);
        end
    end
    
    % We now obtain the routing matrix P by ignoring the non-stateful
    % nodes and calculating by the stochastic complement method the
    % correct transition probabilities, that includes the effects
    % of the non-stateful nodes (e.g., ClassSwitch)
    statefulNodesClasses = [];
    for i=self.getIndexStatefulNodes()
        statefulNodesClasses(end+1:end+K)= ((i-1)*K+1):(i*K);
    end
    
    % Hide the nodes that are not stateful
    rt = dtmc_stochcomp(rtNodes,statefulNodesClasses);
    if nargout >= 4
        rtNodesByClass = cellzeros(K,K,M,M);
        for i=1:M
            for j=1:M
                for r=1:K
                    for s=1:K
                        rtNodesByClass{s,r}(i,j) = rtNodes((i-1)*K+s,(j-1)*K+r);
                    end
                end
            end
        end
    end
    
    if nargout >= 5
        rtNodesByStation = cellzeros(M,M,K,K);
        for i=1:M
            for j=1:M
                for r=1:K
                    for s=1:K
                        rtNodesByStation{i,j}(r,s) = rtNodes((i-1)*K+s,(j-1)*K+r);
                    end
                end
            end
        end
    end
end
