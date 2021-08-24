function isValid = isValid(qn, n, s, options)
% ISVALID = ISVALID(QN, N, S, OPTIONS)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

% n(r): number of jobs at the station in class r
% s(r): jobs of class r that are running

if nargin==3
    options = Solver.defaultOptions;
end
isValid = true;
if isa(qn,'Network')
    qn=qn.getStruct();
end

if isempty(n) & ~isempty(s)
    isValid = false;
    return
end

if iscell(n) %then n is a cell array of states
    ncell = n;
    n = [];
    for isf=1:length(ncell)
        ist = qn.statefulToStation(isf);
        ind = qn.statefulToNode(isf);
        [~, n(ist,:), s(ist,:), ~] = State.toMarginal(qn, ind, ncell{isf}, options);
    end
end

R = qn.nclasses;
K = zeros(1,R);
for ist=1:qn.nstations
    for r=1:R
        K(r) = qn.phases(ist,r);
        if ~isempty(qn.proc) && ~isempty(qn.proc{ist,r}) && any(any(isnan(qn.proc{ist,r}{1}))) && n(ist,r)>0 % if disabled
            isValid = false;
            %            line_error(mfilename,'Chain %d is initialized with an incorrect number of jobs: %f instead of %d.', nc, statejobs_chain, njobs_chain);
            return
        end
    end
    if any(n(ist,:)>qn.classcap(ist,:))
        isValid = false;
        return
    end
end

if nargin > 2 && ~isempty(s)
    for ist=1:qn.nstations
        if qn.nservers(ist)>0
            % if more running jobs than servers
            if sum(s(ist,:)) > qn.nservers(ist)
                switch qn.sched(ist) % don't flag invalid if ps
                    case {SchedStrategy.FCFS,SchedStrategy.SIRO,SchedStrategy.LCFS,SchedStrategy.HOL}
                        isValid = false;
                        return
                end
            end
            % if more running jobs than jobs at the node
            if any(n<s)
                isValid = false;
                return
            end
            % non-idling condition
            if sum(s(ist,:)) ~= min(sum(n(ist,:)), qn.nservers(ist))
                % commented because in ps queues s are in service as well
                % isValid = false;
                %return
            end
        end
    end
end

for nc=1:qn.nchains
    njobs_chain = sum(qn.njobs(find(qn.chains(nc,:))));
    if ~isinf(njobs_chain)
        statejobs_chain = sum(sum(n(:,find(qn.chains(nc,:))),2),1);
        %if ~options.force && abs(1-njobs_chain/statejobs_chain) > options.iter_tol
        if abs(1-njobs_chain/statejobs_chain) > 1e-4
            isValid = false;
            line_error(mfilename,'Chain %d is initialized with an incorrect number of jobs: %f instead of %d.', nc, statejobs_chain, njobs_chain);
            return
        end
        %end
    end
end
end
