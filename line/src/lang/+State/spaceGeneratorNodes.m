function [nodeStateSpace, qn, capacityc] = spaceGeneratorNodes(qn, cutoff, options)
if nargin<3
    options = Solver.defaultOptions;
end
N = qn.njobs';
qn.space = {};
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
        qn.space{isf} = State.fromMarginalBounds(qn, ind, [], capacityc(ind,:), qn.cap(ist), options);
        if isinf(qn.nservers(ist))
            qn.nservers(ist) = sum(capacityc(ind,:));
        end
    elseif qn.isstateful(ind) % generate state space of other stateful nodes that are not stations
        %ist = qn.nodeToStation(ind);
        isf = qn.nodeToStateful(ind);
        switch qn.nodetype(ind)
            case NodeType.Cache
                for r=1:qn.nclasses % restrict state space generation to immediate events
                    if isnan(qn.varsparam{ind}.pref{r})
                        capacityc(ind,r) =  1; %
                    else
                        capacityc(ind,r) =  1; %
                    end
                end
            otherwise
                capacityc(ind,:) =  1; %
        end
        state_var = State.spaceLocalVars(qn, ind);
        state_bufsrv = State.fromMarginalBounds(qn, ind, [], capacityc(ind,:), 1, options);
        qn.space{isf} = State.decorate(state_bufsrv,state_var); % generate all possible states for local variables
    end
end
nodeStateSpace = qn.space;
end