function [SS,SSh,qn] = spaceGenerator(qn, cutoff, options)
% [SS,SSH,QN] = SPACEGENERATOR(QN, CUTOFF)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

% State space generator
% SS: state space
% SSh: hashed state space
% qn: updated qn
N = qn.njobs';
Np = N;

if ~exist('cutoff','var') & any(isinf(Np)) % if has open classes
    line_error(mfilename,'Unspecified cutoff for open classes in state space generator.');
end

if prod(size(cutoff))==1
    cutoff = cutoff * ones(qn.nstations, qn.nclasses);
end

[~, qn, capacityc] = State.spaceGeneratorNodes(qn, cutoff, options);

%%
isOpenClass = isinf(Np);
isClosedClass = ~isOpenClass;
for r=1:qn.nclasses %cut-off open classes to finite capacity
    if isOpenClass(r)
        Np(r) = max(capacityc(:,r)); % if replaced by sum stateMarg_i can exceed capacity
    end
end

nstatefulp = qn.nstateful - sum(qn.nodetype == NodeType.Source); % M without sources
n = pprod(Np);
chainStationPos=[];
%J = zeros(1,Mns*R);
%n = pprod(n,Np);
while n>=0
    % this is rather inefficient since n ignores chains and thus it can
    % generated state spaces such as
    %   J =
    %
    %      0     0     0     0
    %      0     0     0     1
    %      0     0     1     0
    %      0     1     0     0
    %      1     0     0     0
    %      0     0     0     1
    %      0     0     1     0
    %      0     1     0     0
    %      1     0     0     0
    %      0     0     0     2
    %      0     0     1     1
    %      0     0     2     0
    %      0     1     0     1
    %      1     0     0     1
    %      0     1     1     0
    %      1     0     1     0
    %      0     2     0     0
    %      1     1     0     0
    %      2     0     0     0
    % that are then in the need for a call to unique
    if all(isOpenClass) | (Np(isClosedClass) == n(isClosedClass))
        chainStationPos = [chainStationPos; State.spaceClosedMultiCS(nstatefulp,n,qn.chains)];
    end
    n = pprod(n,Np);
end
chainStationPos = unique(chainStationPos,'rows');

netstates = cell(size(chainStationPos,1), qn.nstateful);
for j=1:size(chainStationPos,1)
    for ind=1:qn.nnodes
        if qn.nodetype(ind) == NodeType.Source
            isf = qn.nodeToStateful(ind);
            state_i = State.fromMarginal(qn,ind,[]);
            netstates{j,isf} = State.getHash(qn,ind,state_i);
        elseif qn.isstation(ind)
            isf = qn.nodeToStateful(ind);
            stateMarg_i = chainStationPos(j,(isf-sum(qn.nodetype(1:ind-1) == NodeType.Source)):nstatefulp:end);
            if any(stateMarg_i > capacityc(ind,:))
                netstates{j,isf} = State.getHash(qn,ind,[]);
            else
                state_i = State.fromMarginal(qn,ind,stateMarg_i);
                netstates{j,isf} = State.getHash(qn,ind,state_i);
            end
        elseif qn.isstateful(ind)
            isf = qn.nodeToStateful(ind);
            stateMarg_i = chainStationPos(j,(isf-sum(qn.nodetype(1:ind-1) == NodeType.Source)):nstatefulp:end);
            state_i = qn.space{isf};
            if any(stateMarg_i > capacityc(ind,:))
                netstates{j,isf} = State.getHash(qn,ind,[]);
            elseif qn.nodetype(ind) == NodeType.Cache                
                cacheClasses = union(qn.varsparam{ind}.hitclass, qn.varsparam{ind}.missclass);
                %if sum(stateMarg_i(cacheClasses)) > 1 || sum(stateMarg_i(setdiff(1:qn.nclasses,cacheClasses))) > 0                                        
                if sum(stateMarg_i(1:qn.nclasses)) > 1
                    netstates{j,isf} = State.getHash(qn,ind,[]);
                else
                    state_i = state_i(findrows(state_i(:,1:length(stateMarg_i)),stateMarg_i),:);
                    netstates{j,isf} = State.getHash(qn,ind,state_i);
                end
            else
                state_i = state_i(findrows(state_i(:,1:length(stateMarg_i)),stateMarg_i),:);
                netstates{j,isf} = State.getHash(qn,ind,state_i);
            end
        end
    end
end

ctr = 0;
%SS = sparse([]);
SS = [];
SSh = [];
for j=1:size(chainStationPos,1)
    % for each network state
    v = {netstates{j,:}};
    % cycle over lattice
    vN = cellfun(@length,v)-1;
    n = pprod(vN);
    while n >=0
        u={}; h={};
        skip = false;
        for isf=1:length(n)
            h{isf} = v{isf}(1+n(isf));
            if h{isf} < 0
                skip = true;
                break
            end
            u{isf} = qn.space{isf}(v{isf}(1+n(isf)),:);
        end
        if skip == false
            ctr = ctr + 1; % do not move
            SS(ctr,:)=cell2mat(u);
            SSh(ctr,:)=cell2mat(h);
        end
        n = pprod(n,vN);
    end
end
[SS,IA] = unique(SS,'rows');
SSh = SSh(IA,:);
end
