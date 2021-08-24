function [pi,SSq,arvRates,depRates]=solver_ssa_hashed(qn,options)
% [PI,SSQ,ARVRATES,DEPRATES]=SOLVER_SSA_HASHED(QN,OPTIONS)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

% by default the jobs are all initialized in the first valid state

if ~isfield(options,'seed')
    options.seed = 23000;
end
Solver.resetRandomGeneratorSeed(options.seed+labindex-1);

%% generate local state spaces
nstations = qn.nstations;
nstateful = qn.nstateful;
%init_nserver = qn.nservers; % restore Inf at delay nodes
R = qn.nclasses;
N = qn.njobs';
sync = qn.sync;
csmask = qn.csmask;

cutoff = options.cutoff;
if prod(size(cutoff))==1
    cutoff = cutoff * ones(qn.nstations, qn.nclasses);
end

%%
Np = N';
capacityc = zeros(qn.nnodes, qn.nclasses);
for ind=1:qn.nnodes
    if qn.isstation(ind) % place jobs across stations
        ist = qn.nodeToStation(ind);
        isf = qn.nodeToStateful(ind);
        for r=1:qn.nclasses %cut-off open classes to finite capacity
            c = find(qn.chains(:,r));
            if ~isempty(qn.visits{c}) && qn.visits{c}(ist,r) == 0
                capacityc(ind,r) = 0;
            elseif ~isempty(qn.proc) && ~isempty(qn.proc{ist,r}) && any(any(isnan(qn.proc{ist,r}{1}))) % disabled
                capacityc(ind,r) = 0;
            else
                if isinf(N(r))
                    capacityc(ind,r) =  min(cutoff(ist,r), qn.classcap(ist,r));
                else
                    capacityc(ind,r) =  sum(qn.njobs(qn.chains(c,:)));
                end
            end
        end
        if isinf(qn.nservers(ist))
            qn.nservers(ist) = sum(capacityc(ind,:));
        end
        qn.cap(ist,:) = sum(capacityc(ind,:));
        qn.classcap(ist,:) = capacityc(ind,:);
    end
end
%%
if any(isinf(Np))
    Np(isinf(Np)) = 0;
end

init_state_hashed = ones(1,nstateful); % pick the first state in space{i}

%%
arvRatesSamples = zeros(options.samples,nstateful,R);
depRatesSamples = zeros(options.samples,nstateful,R);
A = length(sync);
samples_collected = 1;
state = init_state_hashed;
stateCell = cell(nstateful,1);
for ind=1:qn.nnodes
    if qn.isstateful(ind)
        isf = qn.nodeToStateful(ind);
        stateCell{isf} = qn.space{isf}(state(isf),:);
    end
end
tranSync = zeros(samples_collected,1);
tranState = zeros(samples_collected,1+length(init_state_hashed));
tranState(1,:) = [0,init_state_hashed];
local = qn.nnodes+1;
last_node_a = 0;
last_node_p = 0;
for act=1:A
    node_a{act} = sync{act}.active{1}.node;
    node_p{act} = sync{act}.passive{1}.node;
    class_a{act} = sync{act}.active{1}.class;
    class_p{act} = sync{act}.passive{1}.class;
    event_a{act} = sync{act}.active{1}.event;
    event_p{act} = sync{act}.passive{1}.event;
    outprob_a{act} = [];
    outprob_p{act} = [];
end

while samples_collected < options.samples
    %samples_collected
    ctr = 1;
    enabled_sync = {}; % row is action label, col1=rate, col2=new state
    enabled_new_state = {};
    enabled_rates = [];
    for act=1:A
        update_cond_a = ((node_a{act} == last_node_a || node_a{act} == last_node_p));
        
        if update_cond_a || isempty(outprob_a{act})
            state_a(act) = state(qn.nodeToStateful(node_a{act}));
            [new_state_a{act}, rate_a{act}, outprob_a{act}, qn] = State.afterEventHashedOrAdd(qn, node_a{act}, state_a(act), event_a{act}, class_a{act});
        end
        
        if all(new_state_a{act}) == -1 % hash not found
            continue
        end
        
        for ia=1:length(new_state_a{act}) % for all possible new states
            if new_state_a{act}(ia) == -1 % hash not found
                continue
            end
            %update_cond_p = ((node_p{act} == last_node_a || node_p{act} == last_node_p)) || isempty(outprob_p{act});
            update_cond = true; %update_cond_a || update_cond_p;
            if rate_a{act}(ia)>0
                if node_p{act} ~= local
                    state_p{act} = state(qn.nodeToStateful(node_p{act}));
                    if node_p{act} == node_a{act} %self-loop
                        if update_cond
                            [new_state_p{act}, ~, outprob_p{act}] = State.afterEventHashedOrAdd(qn, node_p{act}, new_state_a{act}(ia), event_p{act}, class_p{act});
                        end
                    else % departure
                        if update_cond
                            [new_state_p{act}, ~, outprob_p{act}] = State.afterEventHashedOrAdd( qn, node_p{act}, state_p{act}, event_p{act}, class_p{act});
                        end
                    end
                    if new_state_p{act} ~= -1
                        if qn.isstatedep(node_a{act},3)
                            newStateCell = stateCell;
                            newStateCell{qn.nodeToStateful(node_a{act})} = qn.space{qn.nodeToStateful(node_a{act})}(new_state_a{act}(ia),:);
                            newStateCell{qn.nodeToStateful(node_p{act})} = qn.space{qn.nodeToStateful(node_p{act})}(new_state_p{act},:);
                            prob_sync_p{act} = sync{act}.passive{1}.prob(stateCell, newStateCell); %state-dependent
                        else
                            prob_sync_p{act} = sync{act}.passive{1}.prob;
                        end
                    else
                        prob_sync_p{act} = 0;
                    end
                end
                if ~isempty(new_state_a{act}(ia))
                    if node_p{act} == local
                        new_state = state;
                        new_state(qn.nodeToStateful(node_a{act})) = new_state_a{act}(ia);
                        prob_sync_p{act} = 1;
                    elseif new_state_p{act} ~= -1
                        new_state = state;
                        new_state(qn.nodeToStateful(node_a{act})) = new_state_a{act}(ia);
                        new_state(qn.nodeToStateful(node_p{act})) = new_state_p{act};
                    end
                    if ~isnan(rate_a{act})
                        if all(new_state>0)
                            if event_a{act} == EventType.DEP
                                node_a_sf{act} = qn.nodeToStateful(node_a{act});
                                node_p_sf{act} = qn.nodeToStateful(node_p{act});
                                depRatesSamples(samples_collected,node_a_sf{act},class_a{act}) = depRatesSamples(samples_collected,node_a_sf{act},class_a{act}) + outprob_a{act} * outprob_p{act} * rate_a{act}(ia) * prob_sync_p{act};
                                arvRatesSamples(samples_collected,node_p_sf{act},class_p{act}) = arvRatesSamples(samples_collected,node_p_sf{act},class_p{act}) + outprob_a{act} * outprob_p{act} * rate_a{act}(ia) * prob_sync_p{act};
                            end
                            % simulate also self-loops as we need to log them
                            %if any(new_state ~= state)
                            if node_p{act} < local && ~csmask(class_a{act}, class_p{act}) && qn.nodetype(node_p{act})~=NodeType.Source && (rate_a{act}(ia) * prob_sync_p{act} >0)
                                line_error(mfilename,'Fatal error: state-dependent routing at node %d violates class switching mask (node %d -> node %d, class %d -> class %d).', node_a{act}, node_a{act}, node_p{act}, class_a{act}, class_p{act});
                            end
                            enabled_rates(ctr) = rate_a{act}(ia) * prob_sync_p{act};
                            enabled_sync{ctr} = act;
                            enabled_new_state{ctr} = new_state;
                            ctr = ctr + 1;
                            %end
                        end
                    end
                end
            end
        end
    end
    tot_rate = sum(enabled_rates);
    cum_rate = cumsum(enabled_rates) / tot_rate;
    firing_ctr = 1 + max([0,find( rand > cum_rate )]); % select action
    last_node_a = node_a{enabled_sync{firing_ctr}};
    last_node_p = node_p{enabled_sync{firing_ctr}};
    next_state = enabled_new_state{firing_ctr};
    tranState(1+samples_collected, :) = [-(log(rand)/tot_rate), state];
    tranSync(ctr) = enabled_sync{firing_ctr};
    
    samples_collected = samples_collected + 1;
    state = next_state;
    if options.verbose
        if samples_collected == 1e2
            line_printf(sprintf('\b\nSSA samples: %6d',samples_collected));
        elseif options.verbose == 2
            if samples_collected == 0
                line_printf(sprintf('\b\nSSA samples: %6d',samples_collected));
            else
                line_printf(sprintf('\b\b\b\b\b\b\b%6d',samples_collected));
            end
        elseif mod(samples_collected,1e2)==0 || options.verbose == 2
            line_printf(sprintf('\b\b\b\b\b\b\b%6d',samples_collected));
        end
    end
end

%Q
%transient = min([floor(samples_collected/10),1000]); % remove first part of simulation (10% of the samples up to 1000 max)
%transient = 0;
%output = output((transient+1):end,:);
[u,ui,uj] = unique(tranState(:,2:end),'rows');
arvRates = zeros(size(u,1),qn.nstateful,R);
depRates = zeros(size(u,1),qn.nstateful,R);
pi = zeros(1,size(u,1));
for s=1:size(u,1)
    pi(s) = sum(tranState(uj==s,1));
end

for ind=1:qn.nnodes
    if qn.isstateful(ind)
        isf = qn.nodeToStateful(ind);
        if qn.isstation(ind)
            ist = qn.nodeToStation(ind);
            K = qn.phasessz(ist,:);% disabled classes still occupy 1 state element
            Ks = qn.phaseshift(ist,:);
        end
        for s=1:size(u,1)
            for r=1:R
                arvRates(s,isf,r) = arvRatesSamples(ui(s),isf,r); % we just need one sample
                depRates(s,isf,r) = depRatesSamples(ui(s),isf,r); % we just need one sample
            end
            state_i = qn.space{isf}(u(s,isf),:);
            %if isf==2
            %    state_i
            %end
            if qn.isstation(ind)
                [~,nir] = State.toMarginalAggr(qn, ind, state_i, K, Ks);
                ist = qn.nodeToStation(ind);
                SSq(s,((ist-1)*R+1):ist*R) = nir;
            end
        end
    end
end

pi = pi/sum(pi);
%unique(Q,'rows')
%qn.nservers = init_nserver; % restore Inf at delay nodes
end
