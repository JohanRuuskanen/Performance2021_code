function self = link(self, P)
% SELF = LINK(P)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

isReset = false;
if ~isempty(self.qn)
    isReset = true;
    self.resetNetwork; % remove artificial class switch nodes
end
R = self.getNumberOfClasses;
M = self.getNumberOfNodes;

if ~iscell(P) && R>1
    line_error(mfilename,'Multiclass model: the linked routing matrix P must be a cell array, e.g., P = model.initRoutingMatrix; P{1} = Pclass1; P{2} = Pclass2.');
end

isLinearP = true;
if size(P,1) == size(P,2)
    for s=2:R
        for r=1:R
            if nnz(P{r,s})>0
                isLinearP = false;
                break;
            end
        end
    end
    % in this case it is possible that P is linear but just because the
    % routing is state-dependent and therefore some zero entries are
    % actually unspecified
    %cacheNodes = find(cellfun(@(c) isa(c,'Cache'), self.getStatefulNodes));
    for ind=1:M
        switch class(self.nodes{ind})
            case 'Cache'
                % note that since a cache needs to distinguish hits and
                % misses, it needs to do class-switch unless the model is
                % degenerate
                isLinearP = false;
                if self.nodes{ind}.server.hitClass == self.nodes{ind}.server.missClass
                    line_warning(mfilename,'Ambiguous use of hitClass and missClass at cache, it is recommended to use different classes.');
                end
        end
    end
end

for i=self.getDummys
    for r=1:R
        if iscell(P)
            if isLinearP
                P{r}(i,self.getSink) = 1.0;
            else
                P{r,r}(i,self.getSink) = 1.0;
            end
        else
            P(i,self.getSink) = 0.0;
        end
    end
end

% This block is to make sure that P = model.initRoutingMatrix; P{2} writes
% into P{2,2} rather than being interpreted as P{2,1}.
if isLinearP
    Ptmp = P;
    P = cell(R,R);
    for r=1:R
        if iscell(Ptmp)
            P{r,r} = Ptmp{r};
        else
            P{r,r} = Ptmp;
        end
        for s=1:R
            if s~=r
                P{r,s} = 0*Ptmp{r};
            end
        end
    end
end

% assign routing for self-looping jobs
for r=1:R
    if isa(self.classes{r},'SelfLoopingClass')
        for s=1:R
            P{r,s} = 0 * P{r,s};
        end
        P{r,r}(self.classes{r}.reference, self.classes{r}.reference) = 1.0;
    end
end

% link virtual sinks automatically to sink
ispool = cellisa(self.nodes,'Sink');
if sum(ispool) > 1
    line_error(mfilename,'The model can have at most one sink node.');
end

if sum(cellisa(self.nodes,'Source')) > 1
    line_error(mfilename,'The model can have at most one source node.');
end
ispool_nnz = find(ispool)';


if ~iscell(P)
    if R>1
        newP = cell(1,R);
        for r=1:R
            newP{r} = P;
        end
        P = newP;
    else %R==1
        % single class
        for i=ispool_nnz
            P((i-1)*R+1:i*R,:)=0;
        end
        Pmat = P;
        P = cell(R,R);
        for r=1:R
            for s=1:R
                P{r,s} = zeros(M);
                for i=1:M
                    for j=1:M
                        P{r,s}(i,j) = Pmat((i-1)*R+r,(j-1)*R+s);
                    end
                end
            end
        end
    end
end

if numel(P) == R
    % 1 matrix per class
    for r=1:R
        for i=ispool_nnz
            P{r}((i-1)*R+1:i*R,:)=0;
        end
    end
    Pmat = P;
    P = cell(R,R);
    for r=1:R
        P{r,r} = Pmat{r};
        for s=setdiff(1:R,r)
            P{r,s} = zeros(M);
        end
    end
end


isemptyP = false(R,R);
for r=1:R
    for s=1:R
        if isempty(P{r,s})
            isemptyP(r,s)= true;
            P{r,s} = zeros(M);
        else
            for i=ispool_nnz
                P{r,s}(i,:)=0;
            end
        end
    end
end

csMatrix = cell(M,M);
for i=1:M
    for j=1:M
        csMatrix{i,j} = zeros(R);
    end
end

for r=1:R
    for s=1:R
        if ~isemptyP(r,s)
            [I,J] = find(P{r,s});
            for k=1:size(I,1)
                csMatrix{I(k),J(k)}(r,s) = P{r,s}(I(k),J(k));
            end
        end
    end
end


%             for r=1:R
%                 Psum=cellsum({P{r,:}})*ones(M,1);
%                 if min(Psum)<1-1e-4
%                   line_error(mfilename,'Invalid routing probabilities (Node %d departures, switching from class %d).',minpos(Psum),r);
%                 end
%                 if max(Psum)>1+1e-4
%                   line_error(sprintf('Invalid routing probabilities (Node %d departures, switching from class %d).',maxpos(Psum),r));
%                 end
%             end


self.linkedRoutingTable = P;

% As we will now create a CS for each link i->j,
% we now condition on the job going from node i to j
for i=1:M
    for j=1:M
        for r=1:R
            S = sum(csMatrix{i,j}(r,:));
            if S>0
                csMatrix{i,j}(r,:)=csMatrix{i,j}(r,:)/S;
            else
                csMatrix{i,j}(r,r)=1.0;
            end
        end
    end
end

csid = zeros(M);
nodeNames = self.getNodeNames;
for i=1:M
    for j=1:M
        if ~isdiag(csMatrix{i,j})
            self.nodes{end+1} = ClassSwitch(self, sprintf('CS_%s_to_%s',nodeNames{i},nodeNames{j}),csMatrix{i,j});
            csid(i,j) = length(self.nodes);
        end
    end
end

Mplus = length(self.nodes); % number of nodes after addition of cs nodes

% resize matrices
for r=1:R
    for s=1:R
        P{r,s}((M+1):Mplus,(M+1):Mplus)=0;
    end
end

for i=1:M
    for j=1:M
        if csid(i,j)>0
            % re-route
            for r=1:R
                for s=1:R
                    if P{r,s}(i,j)>0
                        P{r,r}(i,csid(i,j)) = P{r,r}(i,csid(i,j)) + P{r,s}(i,j);
                        P{r,s}(i,j) = 0;
                    end
                    P{s,s}(csid(i,j),j) = 1;
                end
            end
        end
    end
end

connected = zeros(Mplus);
for r=1:R
    [I,J,S] = find(P{r,r});
    for k=1:length(I)
        if connected(I(k),J(k)) == 0
            self.addLink(self.nodes{I(k)}, self.nodes{J(k)});
            connected(I(k),J(k)) = 1;
        end
        self.nodes{I(k)}.setProbRouting(self.classes{r}, self.nodes{J(k)}, S(k));
    end
end

self.modifiedRoutingTable = P;

% check if the probability out of any node sums to >1.0
pSum = cellsum(P);
isAboveOne = pSum > 1.0 + Distrib.Zero;
if any(isAboveOne)
    for i=find(isAboveOne)
        if self.nodes{i}.schedStrategy ~= SchedStrategy.FORK
            line_error(mfilename,'The total routing probability for jobs leaving node %s in class %s is greater than 1.0.',self.nodes{i}.name,self.classes{r}.name);
        end
        %        elseif pSum < 1.0 - Distrib.Zero % we cannot check this case as class r may not reach station i, in which case its outgoing routing prob is zero
        %            if self.nodes{i}.schedStrategy ~= SchedStrategy.EXT % if not a sink
        %                line_error(mfilename,'The total routing probability for jobs leaving node %s in class %s is less than 1.0.',self.nodes{i}.name,self.classes{r}.name);
        %            end
    end
end

for i=1:M
    if isa(self.nodes{i},'Place')
        self.nodes{i}.init;
    end
end

if isReset
    self.refreshStruct; % without this exception with linkAndLog
end

end
