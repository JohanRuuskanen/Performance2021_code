function [Q,SS,SSq,Dfilt,arvRates,depRates,qn]=solver_ctmc(qn,options)
% [Q,SS,SSQ,DFILT,ARVRATES,DEPRATES,QN]=SOLVER_CTMC(QN,OPTIONS)
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

%% generate state space
%nnodes = qn.nnodes;
nstateful = qn.nstateful;
nclasses = qn.nclasses;
sync = qn.sync;
csmask = qn.csmask;
if isoctave
    %warning off;
end

[SS,SSh,qnc] = State.spaceGenerator(qn.copy, options.cutoff, options);
qn.space = qnc.space;
if options.verbose
    line_printf('\nCTMC state space size: %d states. ',size(SS,1));
end
if ~isfield(options, 'hide_immediate')
    options.hide_immediate = true;
end

%size(SS)
%%
Q = sparse(eye(size(SSh,1))); % the diagonal elements will be removed later
A = length(sync);
Dfilt = cell(1,A);
for a=1:A
    Dfilt{a} = 0*Q;
end
local = qn.nnodes+1; % passive action
arvRates = zeros(size(SSh,1),nstateful,nclasses);
depRates = zeros(size(SSh,1),nstateful,nclasses);
SSq = zeros(size(SSh));

% for all synchronizations
for a=1:A
    stateCell = cell(nstateful,1);
    for s=1:size(SSh,1)
        state = SSh(s,:);
        % update state cell array and SSq
        for ind = 1:qn.nnodes
            if qn.isstateful(ind)
                isf = qn.nodeToStateful(ind);
                stateCell{isf} = qn.space{isf}(state(isf),:);
                if qn.isstation(ind)
                    ist = qn.nodeToStation(ind);
                    [~,nir] = State.toMarginal(qn,ind,stateCell{isf});
                    SSq(s,((ist-1)*nclasses+1):ist*nclasses) = nir;
                end
            end
        end
    end
end

%% for all synchronizations
for a=1:A
    stateCell = cell(nstateful,1);
    %qn.sync{a}.active{1}.print
    for s=1:size(SSh,1)
        %[a,s]
        state = SSh(s,:);
        % update state cell array and SSq
        for ind = 1:qn.nnodes
            if qn.isstateful(ind)
                isf = qn.nodeToStateful(ind);
                stateCell{isf} = qn.space{isf}(state(isf),:);
                if qn.isstation(ind)
                    ist = qn.nodeToStation(ind);
                    [~,nir] = State.toMarginal(qn,ind,stateCell{isf});
                end
            end
        end
        node_a = sync{a}.active{1}.node;
        state_a = state(qn.nodeToStateful(node_a));
        class_a = sync{a}.active{1}.class;
        event_a = sync{a}.active{1}.event;
        [new_state_a, rate_a] = State.afterEventHashed( qn, node_a, state_a, event_a, class_a);        
        %% debugging block
%        if true%options.verbose == 2
%            line_printf('---\n');
%            sync{a}.active{1}.print,
%        end
        %%
        if new_state_a == -1 % hash not found
            continue
        end
        for ia=1:length(new_state_a)
            if rate_a(ia)>0
                node_p = sync{a}.passive{1}.node;
                if node_p ~= local
                    state_p = state(qn.nodeToStateful(node_p));
                    class_p = sync{a}.passive{1}.class;
                    event_p = sync{a}.passive{1}.event;
                    %prob_sync_p = sync{a}.passive{1}.prob(state_a, state_p)
                    %if prob_sync_p > 0
                    %% debugging block
                    if options.verbose == 2
                        line_printf('---\n');
                        sync{a}.active{1}.print,
                        sync{a}.passive{1}.print
                    end
                    %%
                    if node_p == node_a %self-loop
                        [new_state_p, ~, outprob_p] = State.afterEventHashed( qn, node_p, new_state_a(ia), event_p, class_p);
                    else % departure
                        [new_state_p, ~, outprob_p] = State.afterEventHashed( qn, node_p, state_p, event_p, class_p);
                    end
                    for ip=1:size(new_state_p,1)
                        if node_p ~= local
                            if new_state_p ~= -1
                                if qn.isstatedep(node_a,3)
                                    newStateCell = stateCell;
                                    newStateCell{qn.nodeToStateful(node_a)} = qn.space{qn.nodeToStateful(node_a)}(new_state_a(ia),:);
                                    newStateCell{qn.nodeToStateful(node_p)} = qn.space{qn.nodeToStateful(node_p)}(new_state_p(ip),:);
                                    prob_sync_p = sync{a}.passive{1}.prob(stateCell, newStateCell) * outprob_p(ip); %state-dependent
                                else
                                    prob_sync_p = sync{a}.passive{1}.prob * outprob_p(ip);
                                end
                            else
                                prob_sync_p = 0;
                            end
                        end
                        if ~isempty(new_state_a(ia))
                            if node_p == local % local action
                                new_state = state;
                                new_state(qn.nodeToStateful(node_a)) = new_state_a(ia);
                                prob_sync_p = outprob_p(ip);
                            elseif ~isempty(new_state_p)
                                new_state = state;
                                new_state(qn.nodeToStateful(node_a)) = new_state_a(ia);
                                new_state(qn.nodeToStateful(node_p)) = new_state_p(ip);
                            end
                            ns = matchrow(SSh, new_state);
                            if ns>0
                                if ~isnan(rate_a)
                                    if node_p < local && ~csmask(class_a, class_p) && rate_a(ia) * prob_sync_p >0 && (qn.nodetype(node_p)~=NodeType.Source)
                                        line_error(mfilename,'Error: state-dependent routing at node %d (%s) violates the class switching mask (node %d -> node %d, class %d -> class %d).', node_a, qn.nodenames{node_a}, node_a, node_p, class_a, class_p);
                                    end
                                    if size(Dfilt{a}) >= [s,ns] % check needed as D{a} is a sparse matrix
                                        Dfilt{a}(s,ns) = Dfilt{a}(s,ns) + rate_a(ia) * prob_sync_p;
                                    else
                                        Dfilt{a}(s,ns) = rate_a(ia) * prob_sync_p;
                                    end
                                end
                            end
                        end
                    end
                else % node_p == local
                    if ~isempty(new_state_a(ia))
                        new_state = state;
                        new_state(qn.nodeToStateful(node_a)) = new_state_a(ia);
                        prob_sync_p = 1;
                        ns = matchrow(SSh, new_state);
                        if ns>0
                            if ~isnan(rate_a)
                                if size(Dfilt{a}) >= [s,ns] % needed for sparse matrix
                                    Dfilt{a}(s,ns) = Dfilt{a}(s,ns) + rate_a(ia) * prob_sync_p;
                                else
                                    Dfilt{a}(s,ns) = rate_a(ia) * prob_sync_p;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%%
for a=1:A
    Q = Q + Dfilt{a};
    % active
    node_a = sync{a}.active{1}.node;
    class_a = sync{a}.active{1}.class;
    event_a = sync{a}.active{1}.event;
    % passive
    node_p = sync{a}.passive{1}.node;
    class_p = sync{a}.passive{1}.class;
    if event_a == EventType.DEP
        node_a_sf = qn.nodeToStateful(node_a);
        node_p_sf = qn.nodeToStateful(node_p);
        for s=1:size(SSh,1)
            depRates(s,node_a_sf,class_a) = depRates(s,node_a_sf,class_a) + sum(Dfilt{a}(s,:));
            arvRates(s,node_p_sf,class_p) = arvRates(s,node_p_sf,class_p) + sum(Dfilt{a}(s,:));
        end
    end
end
zero_row = find(sum(Q,2)==0);
zero_col = find(sum(Q,1)==0);

%%
% in case the last column of Q represent a state for a transient class, it
% is possible that no transitions go back to it, although it is valid for
% the system to be initialized in that state. So we need to fill-in the
% zeros at the end.
Q(:,end+1:end+(size(Q,1)-size(Q,2)))=0;
Q(zero_row,zero_row) = -eye(length(zero_row)); % can this be replaced by []?
Q(zero_col,zero_col) = -eye(length(zero_col));
for a=1:A
    Dfilt{a}(:,end+1:end+(size(Dfilt{a},1)-size(Dfilt{a},2)))=0;
end

if options.verbose == 2
    SolverCTMC.printInfGen(Q,SS);
end
Q = ctmc_makeinfgen(Q);

%SolverCTMC.printInfGen(Q,SS)
%% now remove immediate transitions
% we first determine states in stateful nodes where there is an immediate
% job in the node
if options.hide_immediate % if want to remove immediate transitions
    statefuls = find(qn.isstateful-qn.isstation);
    imm=[];
    for st = statefuls(:)'
        imm_st = find(sum(qn.space{st}(:,1:nclasses),2)>0);
        imm = [imm; find(arrayfun(@(a) any(a==imm_st),SSh(:,st)))];
    end
    imm = unique(imm);
    nonimm = setdiff(1:size(Q,1),imm);
    SS(imm,:) = [];
    SSq(imm,:) = [];
    arvRates(imm,:,:) = [];
    depRates(imm,:,:) = [];
    [Q,~,~,~,~,T] = ctmc_stochcomp(Q, nonimm);
    for a=1:A
        Dfilt{a} = Dfilt{a}(nonimm,nonimm)+T;
    end
end
%SolverCTMC.printInfGen(Q,SS)
%%
Q = ctmc_makeinfgen(Q);
end
