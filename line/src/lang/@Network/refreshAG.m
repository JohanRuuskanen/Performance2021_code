function refreshAG(self, cutoffs)
% REFRESHAG()
% Export network as a collection of interacting stochastic agents


% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

% We now generate the global state space
% parses all but the service processes
qn = self.getStruct;
nnodes = qn.nnodes;
nstateful = qn.nstateful;
S = length(qn.sync);
%%
if nargin>1 && length(cutoffs(:))==1
    defaultCutoff = cutoffs;
elseif nargin==1 % cutoff not specified, use 2 in each dimension
    defaultCutoff = 2;
end
if nargin==1 || (nargin>1 && length(cutoffs(:))==1)
    cutoffs = defaultCutoff * ones(qn.nnodes,qn.nclasses);
end
qn = self.getStruct;
space = State.spaceGeneratorNodes(qn, cutoffs);
%%
active = 1; % column of active actions
passive = 2; % columns of passive actions
eps = nstateful+1; % row of local actions
sync = qn.sync;

RL = cell(S+1,nstateful); % RL{a,active}(s1,s2)=rate from state s1 to s2 for action a in active agent, RL{a,passive}(s1p,s2p)=accept probability of passive

for j=1:nnodes
    if qn.isstateful(j)
        jsf = qn.nodeToStateful(j);
        RL{end,jsf}=zeros(size(space{jsf},1));
    end
end

issim = false;
for s=1:S
    node_a = sync{s}.active{1}.node;
    class_a = sync{s}.active{1}.class;
    event_a = sync{s}.active{1}.event;
    node_p = sync{s}.passive{1}.node;
    class_p = sync{s}.passive{1}.class;
    event_p = sync{s}.passive{1}.event; if isempty(event_p), event_p=0; end
    prob_p = sync{s}.passive{1}.prob;
    
    if node_p ~= eps
        if prob_p > 0
            %a = matchrow(actionTable,[node_a, class_a, node_p, class_p, event_a, event_p]);
            % determine active-passive roles
            AP(s,active) = qn.nodeToStateful(node_a); % active
            AP(s,passive) = qn.nodeToStateful(node_p);
            % determine synchronization rates
            RL{s,active} = zeros(size(space{qn.nodeToStateful(node_a)},1));
            RL{s,passive} = zeros(size(space{qn.nodeToStateful(node_p)},1));
            for sa=1:size(space{qn.nodeToStateful(node_a)},1)
                state_a = space{qn.nodeToStateful(node_a)}(sa,:);
                [new_state_a, rate_a] = State.afterEvent(qn, node_a, state_a, event_a, class_a, issim);
                if ~isempty(new_state_a)
                    nsa = matchrow(space{qn.nodeToStateful(node_a)},new_state_a);
                    if sa>0 && nsa >= 0
                        RL{s,active}(sa,nsa)=rate_a;
                        
                        node_p = sync{s}.passive{1}.node;
                        state_p = space{qn.nodeToStateful(node_p)};
                        if node_p == node_a %self-loop
                            [new_state_p, ~] = State.afterEvent(qn, node_p, new_state_a, event_p, class_p, issim);
                        else
                            [new_state_p, ~] = State.afterEvent(qn, node_p, state_p, event_p, class_p, issim);
                        end
                        for j = 1:size(state_p,1)
                            try
                                sp = matchrow(space{qn.nodeToStateful(node_p)},state_p(j,:));
                                nsp = matchrow(space{qn.nodeToStateful(node_p)},new_state_p(j,:));
                                if sp>0 && nsp >= 0
                                    RL{s,passive}(sp,nsp)=prob_p;
                                end
                            end
                        end
                    end
                end
            end
        end
    else % local transitions
        for sa=1:size(space{qn.nodeToStateful(node_a)},1)
            state_a = space{qn.nodeToStateful(node_a)}(sa,:);
            [new_state_a, rate_a] = State.afterEvent(qn, node_a, state_a, event_a, class_a, issim);
            nsa = matchrow(space{qn.nodeToStateful(node_a)},new_state_a);
            if sa > 0 && nsa > 0
                RL{end,qn.nodeToStateful(node_a)}(sa,nsa)=RL{end,qn.nodeToStateful(node_a)}(sa,nsa)+rate_a;
            end
        end
    end
end

% remove inactive synchornizations
inactiveSync = [];
for a=1:size(RL,1)-1 % don't delete local events row
    if ~any(RL{a,1}(:))
        inactiveSync = [inactiveSync a];
    end
end
RL(inactiveSync,:)=[];
AP(inactiveSync,:)=[];
self.ag = struct();
self.ag.space = space;
self.ag.actions = RL;
self.ag.role = AP;
%[x, pi, Q, stats, xls, xus] = autocat(ag.rate, ag.role, 'linprog', 'lpr')
end
