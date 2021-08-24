function space = fromMarginal(qn, ind, n, options)
% SPACE = FROMMARGINAL(QN, IND, N, OPTIONS)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

if ~exist('options','var')
    options.force = false;
end
if isa(qn,'Network')
    qn=qn.getStruct();
end
% generate states such that the marginal queue-lengths are as in vector n
%  n(r): number of jobs at the station in class r
R = qn.nclasses;
S = qn.nservers;
state = [];
space = [];

% ind: node index
ist = qn.nodeToStation(ind);
isf = qn.nodeToStateful(ind);

if qn.isstateful(ind) && ~qn.isstation(ind)
    for r=1:R
        init = State.spaceClosedSingle(1,n(r));
        state = State.decorate(state,init);
    end
    space = State.decorate(space,state);
    return
end

K = zeros(1,R);
for r=1:R
    if isempty(qn.proc{ist,r})
        K(r) = 0;
    else
        K(r) = length(qn.proc{ist,r}{1});
    end
end
if (qn.schedid(ist) ~= SchedStrategy.ID_EXT) && any(n>qn.classcap(ist,:))
    return
end

% generate local-state space
switch qn.nodetype(ind)
    case {NodeType.Queue, NodeType.Delay, NodeType.Source}
        switch qn.sched(ist)
            case SchedStrategy.EXT
                for r=1:R
                    if ~isempty(qn.proc) && ~isempty(qn.proc{ist,r}) && any(any(isnan(qn.proc{ist,r}{1}))) % disabled
                        init = 0*ones(1,K(r));
                    else
                        init = State.spaceClosedSingle(K(r),1);
                    end
                    state = State.decorate(state,init);
                end
                space = State.decorate(space,state);
                space = [Inf*ones(size(space,1),1),space];
            case {SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS}
                % in these policies we only track the jobs in the servers
                for r=1:R
                    init = State.spaceClosedSingle(K(r),n(r));
                    state = State.decorate(state,init);
                end
                space = State.decorate(space,state);
            case {SchedStrategy.SIRO, SchedStrategy.LEPT, SchedStrategy.SEPT}
                % in these policies we track an un-ordered buffer and
                % the jobs in the servers
                % build list of job classes in the node, with repetition
                if sum(n) <= S(ist)
                    for r=1:R
                        init = State.spaceClosedSingle(K(r),n(r));
                        state = State.decorate(state,init);
                    end
                    space = State.decorate(space,[zeros(size(state,1),R),state]);
                else
                    si = multichoosecon(n,S(ist)); % jobs of class r that are running
                    mi_buf = repmat(n,size(si,1),1) - si; % jobs of class r in buffer
                    for k=1:size(si,1)
                        % determine number of classes r jobs running in phase j
                        kstate=[];
                        for r=1:R
                            init = State.spaceClosedSingle(K(r),si(k,r));
                            kstate = State.decorate(kstate,init);
                        end
                        state = [repmat(mi_buf(k,:),size(kstate,1),1), kstate];
                        space = [space; state];
                    end
                end
            case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.LCFS}
                sizeEstimator = multinomialln(n) - gammaln(sum(n)) + gammaln(1+qn.cap(ist));
                sizeEstimator = round(sizeEstimator/log(10));
                if sizeEstimator > 2
                    if ~isfield(options,'force') || options.force == false
                        %line_warning(mfilename,sprintf('fromMarginal(): Marginal state space size is very large: 1e%d states. Set options.force=true to bypass this control.\n',sizeEstimator));
                    end
                end
                
                if sum(n) == 0
                    space = zeros(1,1+sum(K)); % unclear if this should 1+sum(K), was sum(K) but State.fromMarginalAndStarted uses 1+sum(K) so was changed here as well
                    return
                end
                % in these policies we track an ordered buffer and
                % the jobs in the servers
                
                % build list of job classes in the node, with repetition
                vi = [];
                for r=1:R
                    if n(r)>0
                        vi=[vi, r*ones(1,n(r))];
                    end
                end
                
                % gen permutation of their positions in the fcfs buffer
                mi = uniqueperms(vi);
                if isempty(mi)
                    mi_buf = zeros(1,max(0,sum(n)-S(ist)));
                    state = zeros(1,R);
                    state = State.decorate(state,[mi_buf,state]);
                else                    
                    mi = mi(:,(end-min(sum(n),qn.cap(ist))+1):end); % n(r) may count more than once elements within the same chain
                    mi = unique(mi,'rows');
                    % mi_buf: class of job in buffer position i (0=empty)
                    mi_buf = [zeros(size(mi,1),min(sum(n),qn.cap(ist))-S(ist)-size(mi(:,1:end-S(ist)),2)), mi(:,1:end-S(ist))];
                    if isempty(mi_buf)
                        mi_buf = zeros(size(mi,1),1);
                    end
                    % mi_srv: class of job running in server i
                    mi_srv = mi(:,max(size(mi,2)-S(ist)+1,1):end);
                    % si: number of class r jobs that are running
                    si =[];
                    for k=1:size(mi_srv,1)
                        si(k,1:R) = hist(mi_srv(k,:),1:R);
                    end
                    %si = unique(si,'rows');
                    for k=1:size(si,1)
                        % determine number of class r jobs running in phase
                        % j in server state mi_srv(kjs,:) and build
                        % state
                        kstate=[];
                        for r=1:R
                            kstate = State.decorate(kstate,State.spaceClosedSingle(K(r),si(k,r)));
                        end
                        state = [state; repmat(mi_buf(k,:),size(kstate,1),1), kstate];
                    end
                end
                space = state;
            case {SchedStrategy.SJF, SchedStrategy.LJF}
                % in these policies the state space includes continuous
                % random variables for the service times
                line_error(mfilename,'The scheduling policy does not admit a discrete state space.\n');
        end
        switch qn.routing(ind)
            case RoutingStrategy.ID_RRB
                space = State.decorate(space, qn.varsparam{ind}.outlinks(:));
        end
    case NodeType.Cache
        switch qn.sched(ist)
            case SchedStrategy.INF
                % in this policies we only track the jobs in the servers
                for r=1:R
                    init = State.spaceClosedSingle(K(r),n(r));
                    state = State.decorate(state,init);
                end
                space = State.decorate(space,state);
        end
        switch qn.routing(ind)
            case RoutingStrategy.ID_RRB
                space = State.decorate(space, qn.varsparam{ind}.outlinks(:));
        end
end
space = unique(space,'rows'); % do not comment, required to sort empty state as first
space = space(end:-1:1,:); % so that states with jobs in phase 1 comes earlier
end
