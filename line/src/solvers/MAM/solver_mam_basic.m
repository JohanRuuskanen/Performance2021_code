function [QN,UN,RN,TN,CN,XN] = solver_mam_basic(qn, options, config)
% [Q,U,R,T,C,X] = SOLVER_MAM(QN, PH, OPTIONS)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

global BuToolsVerbose;
global BuToolsCheckInput;
global BuToolsCheckPrecision;

PH = qn.proc;
%% generate local state spaces
I = qn.nnodes;
M = qn.nstations;
K = qn.nclasses;
C = qn.nchains;
N = qn.njobs';
V = cellsum(qn.visits);

QN = zeros(M,K);
UN = zeros(M,K);
RN = zeros(M,K);
TN = zeros(M,K);
CN = zeros(1,K);
XN = zeros(1,K);

lambda = zeros(1,K);
lambdas_inchain = {};
for c=1:C
    inchain = find(qn.chains(c,:));
    lambdas_inchain{c} = qn.rates(qn.refstat(inchain(1)),inchain);
    %lambdas_inchain{c} = lambdas_inchain{c}(isfinite(lambdas_inchain{c}));
    lambda(inchain) = sum(lambdas_inchain{c}(isfinite(lambdas_inchain{c})));
end

if qn.isopen()
    % open queueing system (one node is the external world)
    BuToolsVerbose = false;
    BuToolsCheckInput = true;
    BuToolsCheckPrecision = 1e-12;
    pie = {};
    D0 = {};
    % first build the joint arrival process
    for ist=1:M
        switch qn.sched(ist)
            case SchedStrategy.EXT
                % assemble a MMAP for the arrival process from all classes
                for k=1:K
                    if isnan(PH{ist,k}{1})
                        PH{ist,k} = map_exponential(Inf); % no arrivals from this class
                    end
                end
                chainArrivalAtSource = cell(1,C);
                for c=1:C %for each chain
                    inchain = find(qn.chains(c,:))';
                    k = inchain(1);
                    chainArrivalAtSource{c} = {PH{ist,k}{1},PH{ist,k}{2},PH{ist,k}{2}};
                    for ki=2:length(inchain)
                        k = inchain(ki);
                        if isnan(PH{ist,k}{1})
                            PH{ist,k} = map_exponential(Inf); % no arrivals from this class
                        end
                        chainArrivalAtSource{c} = mmap_super_safe({chainArrivalAtSource{c},{PH{ist,k}{1},PH{ist,k}{2},PH{ist,k}{2}}}, config.space_max, 'default');
                    end
%                     if c == 1
%                         aggrArrivalAtSource = mmap_super_safe({chainArrivalAtSource{1}, mmap_exponential(0,1)}, config.space_max, 'default');
%                         aggrArrivalAtSource = {aggrArrivalAtSource{1} aggrArrivalAtSource{2} aggrArrivalAtSource{2}};
%                         aggrArrivalAtSource = mmap_scale(aggrArrivalAtSource, 1/ map_lambda(chainArrivalAtSource{c}));
%                     else
%                         aggrArrivalAtSource = mmap_super_safe({aggrArrivalAtSource, chainArrivalAtSource{c}},config.space_max, 'default');
%                     end
                    inchain = find(qn.chains(c,:));
                    TN(ist,inchain) = lambdas_inchain{c};
                    %TN(ist,isnan(lambdas_inchain{c}))=0;
                end
            case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.PS}
                for k=1:K
                    % divide service time by number of servers and put
                    % later a surrogate delay server in tandem to compensate
                    PH{ist,k} = map_scale(PH{ist,k}, map_mean(PH{ist,k})/qn.nservers(ist));
                    pie{ist,k} = map_pie(PH{ist,k});
                    D0{ist,k} = PH{ist,k}{1};
                end
        end
    end % i
    
    % at the first iteration, propagate the arrivals with the same
    for ind=1:I
        if qn.isstation(ind)
            ist = qn.nodeToStation(ind);
            switch qn.sched(ist)
                case SchedStrategy.INF
                    for k=1:K
                        TN(ist,k) = lambda(k)*V(ist,k);
                        UN(ist,k) = map_mean(PH{ist,k})*TN(ist,k);
                        QN(ist,k) = TN(ist,k).*map_mean(PH{ist,k})*V(ist,k);
                        RN(ist,k) = QN(ist,k)/TN(ist,k);
                    end
                case SchedStrategy.PS
                    for k=1:K
                        TN(ist,k) = lambda(k)*V(ist,k);
                        UN(ist,k) = map_mean(PH{ist,k})*TN(ist,k);
                    end
                    sum(UN(ist,:))
                    for k=1:K
                        QN(ist,k) = UN(ist,k)/(1-sum(UN(ist,:)));
                        RN(ist,k) = QN(ist,k)/TN(ist,k);
                    end
                case {SchedStrategy.FCFS, SchedStrategy.HOL}
                    chainArrivalAtNode = cell(1,C);
                    Qret = {};
                    rates = cell(M,C);
                    for c=1:C %for each chain
                        rates{ist,c} = V(ist,:) .* map_lambda(chainArrivalAtSource{c});
                        inchain = find(qn.chains(c,:))';
                        chainArrivalAtNode{c} = mmap_mark(chainArrivalAtSource{c}, rates{ist,c}(inchain) / sum(rates{ist,c}(inchain)));
                        chainArrivalAtNode{c} = mmap_scale(chainArrivalAtNode{c}, 1./rates{ist,c});
                        if c == 1
                            aggrArrivalAtNode = mmap_super_safe({chainArrivalAtNode{c}, mmap_exponential(0,1)}, config.space_max, 'default');                            
                            aggrArrivalAtNode = {aggrArrivalAtNode{1} aggrArrivalAtNode{2} aggrArrivalAtNode{2}};
                            aggrArrivalAtNode = mmap_scale(aggrArrivalAtNode, 1/ map_lambda(chainArrivalAtNode{c}));                            
                        else
                            aggrArrivalAtNode = mmap_super_safe({aggrArrivalAtNode, chainArrivalAtNode{c}}, config.space_max, 'default');
                        end
                    end
                    if strcmp(qn.sched(ist),SchedStrategy.HOL) && any(qn.classprio ~= qn.classprio(1)) % if priorities are not identical
                        [uK,iK] = unique(qn.classprio);
                        if length(uK) == length(qn.classprio) % if all priorities are different
                            [Qret{iK}] = MMAPPH1NPPR({aggrArrivalAtNode{[1;2+iK]}}, {pie{ist,iK}}, {D0{ist,iK}}, 'ncMoms', 1);
                        else
                            line_error(mfilename,'Solver MAM requires either identical priorities or all distinct priorities');
                        end
                    else
                        [Qret{1:K}] = MMAPPH1FCFS({aggrArrivalAtNode{[1,3:end]}}, {pie{ist,:}}, {D0{ist,:}}, 'ncMoms', 1);
                    end
                    QN(ist,:) = cell2mat(Qret);
                    for k=1:K
                        c = find(qn.chains(:,k));
                        TN(ist,k) = rates{ist,c}(k);
                        UN(ist,k) = TN(ist,k) * map_mean(PH{ist,k});
                        % add number of jobs at the surrogate delay server
                        QN(ist,k) = QN(ist,k) + TN(ist,k)*(map_mean(PH{ist,k})*qn.nservers(ist)) * (qn.nservers(ist)-1)/qn.nservers(ist);
                        RN(ist,k) = QN(ist,k) ./ TN(ist,k);
                    end
            end
        else % not a station
            switch qn.nodetype(ind)
                case NodeType.Fork
                    
            end
        end
    end
    CN = sum(RN,1);
else
    line_warning(mfilename,'This model is not supported by SolverMAM yet. Returning with no result.');
end

end
