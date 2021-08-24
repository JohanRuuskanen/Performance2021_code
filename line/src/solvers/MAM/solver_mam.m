function [QN,UN,RN,TN,CN,XN] = solver_mam(qn, options, config)
%[Q,U,R,T,C,X] = SOLVER_MAM(QN, PH, OPTIONS)

%Copyright (c) 2012-2020, Imperial College London
%All rights reserved.
global SCVmam
global BuToolsVerbose;
global BuToolsCheckInput;
global BuToolsCheckPrecision;

PH = qn.proc;
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
for c=1:C
    inchain = find(qn.chains(c,:));
    lambdas_inchain = qn.rates(qn.refstat(inchain(1)),inchain);
    lambdas_inchain = lambdas_inchain(isfinite(lambdas_inchain));
    lambda(inchain) = sum(lambdas_inchain);
end

chain = zeros(1,K);
for k=1:K
    chain(k) = find(qn.chains(:,k));
end

if qn.isopen()
    %    open queueing system (one node is the external world)
    BuToolsVerbose = false;
    BuToolsCheckInput = false;
    BuToolsCheckPrecision = 1e-14;
    
    pie = {};
    D0 = {};
    for ist=1:M
        switch qn.sched(ist)
            case SchedStrategy.EXT
                TN(ist,:) = qn.rates(ist,:);
                TN(ist,isnan(TN(ist,:)))=0;
            case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.PS}
                for k=1:K
                    %                    divide service time by number of servers and put
                    %                    later a surrogate delay server in tandem to compensate
                    PH{ist,k} = map_scale(PH{ist,k}, map_mean(PH{ist,k})/qn.nservers(ist));
                    pie{ist,k} = map_pie(PH{ist,k});
                    D0{ist,k} = PH{ist,k}{1};
                end
        end
    end
    
    it_max = options.iter_max;
    for it=1:it_max
        %it
        %        now estimate arrival processes
        if it == 1
            %            initially form departure processes using scaled service
            DEP = PH;
            for ind=1:M
                for r=1:K
                    ist = qn.nodeToStation(ind);
                    DEP{ind,r} = map_scale(PH{ist,r}, 1 / (lambda(r) * V(ind,r)) );
                end
            end
        end
        
        
        ARV = solver_mam_estflows(qn, DEP, config);
        
        QN_1 = QN;
        
        for ist=1:M
            ind = qn.stationToNode(ist);
            switch qn.nodetype(ind)
                case NodeType.Queue
                    if length(ARV{ind}{1}) > config.space_max
                        line_printf('\nArrival process at node %d is now at %d states. Compressing.',ind,length(ARV{ind}{1}));
                        ARV{ind} = mmap_compress(ARV{ind});
                    end                    
                    [Qret{1:K}, ncDistr] = MMAPPH1FCFS({ARV{ind}{[1,3:end]}}, {pie{ist,:}}, {D0{ist,:}}, 'ncMoms', 1, 'ncDistr',2);
                    for k=1:K
                        QN(ist,k) = sum(Qret{k});
                    end
                    TN(ist,:) = mmap_lambda(ARV{ind});
            end
            for k=1:K
                UN(ist,k) = TN(ist,k) * map_mean(PH{ist,k});
                %add number of jobs at the surrogate delay server
                QN(ist,k) = QN(ist,k) + TN(ist,k)*(map_mean(PH{ist,k})*qn.nservers(ist)) * (qn.nservers(ist)-1)/qn.nservers(ist);
                RN(ist,k) = QN(ist,k) ./ TN(ist,k);
            end
        end
        
        if it >=3 && max(abs(QN(:)-QN_1(:))./QN_1(:)) < options.iter_tol
            break;
        end
        
        for ist=1:M
            ind = qn.stationToNode(ist);
            switch qn.nodetype(ind)
                case NodeType.Queue
                    [Ret{1:2*K}] = MMAPPH1FCFS({ARV{ind}{[1,3:end]}}, {pie{ist,:}}, {D0{ist,:}}, 'stDistrPH');
                    for r=1:K
                        
                        %obtain response time distribution for class r
                        %alpha = Ret{(r-1)*2+1}; T0 = Ret{(r-1)*2+2};
                        %RD = {T0,(-T0)*ones(length(alpha),1)*alpha(:)'}; %PH to MAP
                        %tRD = sum(RD{2},2);
                        %pieRD = map_pie(RD);
                        
                        %define a ph for the arrival process of class r
                        A = mmap_hide(ARV{ind},setdiff(1:K,r));
                        tA = sum(A{2},2);
                        pieA = map_pie(A);
                        
                        %define a ph for the service process of class r
                        S = PH{ist,r};
                        pieS = map_pie(S);
                        tS = sum(S{2},2);
                        
                        %with probability rho the output is the
                        %inter-arrival time+service time, with probability 1-rho is the
                        %response time
                        rho = sum(UN(ist,:));
                        AQ = sum(QN(ist,:)); % queue seen upon arrival with PASTA
                        Afull = AQ/rho; % from AQ = (1-rho)*0 + rho*AQfull
                        pfullonarrival = 1 - (Afull)/(1+Afull);
                        
                        A = mmap_scale(A,map_mean(A)-map_mean(S));
                        
                        zAS = 0*tA*pieS;
                        zSA = 0*tS*pieA;
                        zA = 0*A{2};
                        
                        % this is only for single-class
                        DEP0ir = [S{1}, sum(rho)*tS*pieA;
                            zAS, A{1}];
                        
                        DEP1ir = [(1-sum(rho))*S{2}, zSA;
                            tA*pieS, zA];
                        
                        DEP{ind,r} = map_normalize({DEP0ir,DEP1ir});
                        
                        DEP{ind,r} = map_scale(DEP{ind,r}, 1 / (lambda(r) * V(ind,r)) );
                        SCVd(ind,r) = map_scv(DEP{ind,r});
                        IDCd(ind,r) = map_idc(DEP{ind,r});
                    end
            end
        end
        
        if it>1
            SCVmam=SCVd(2:end);
            %SCVmam
            %IDCd
            %QN
        end
    end
    if options.verbose
        line_printf('\nMAM parametric decomposition completed in %d iterations.',it);
    end
else
    line_warning(mfilename,'This model is not supported by SolverMAM yet. Returning with no result.');
end

end
