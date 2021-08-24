function [XN,UN,QN,RN,TN,CN,qnc]=solver_ssa_analysis_spmd(laboptions, qn, qnc, PH)
% [XN,UN,QN,RN,TN,CN]=SOLVER_SSA_ANALYSIS_SPMD(LABOPTIONS, QN, QNC, PH)

M = qn.nstations;    %number of stations
K = qn.nclasses;    %number of classes

S = qn.nservers;
NK = qn.njobs';  % initial population per class
sched = qn.sched;
spmd
    laboptions.samples = ceil(laboptions.samples / numlabs);
    laboptions.verbose = false;
    switch laboptions.method
        case {'para','parallel'}
            [probSysState,SSq,arvRates,depRates] = solver_ssa(qnc, laboptions);
            qn.space = qnc.space;
        case {'para.hash','parallel.hash'}
            [probSysState,SSq,arvRates,depRates] = solver_ssa_hashed(qnc, laboptions);
            qn.space = qnc.space;
            
    end
    XN = NaN*zeros(1,K);
    UN = NaN*zeros(M,K);
    QN = NaN*zeros(M,K);
    RN = NaN*zeros(M,K);
    TN = NaN*zeros(M,K);
    CN = NaN*zeros(1,K);
    for k=1:K
        refsf = qn.stationToStateful(qn.refstat(k));
        XN(k) = probSysState*arvRates(:,refsf,k);
        for i=1:M
            isf = qn.stationToStateful(i);
            TN(i,k) = probSysState*depRates(:,isf,k);
            QN(i,k) = probSysState*SSq(:,(i-1)*K+k);
            switch qn.sched(i)
                case SchedStrategy.INF
                    UN(i,k) = QN(i,k);
                otherwise
                    % we use Little's law, otherwise there are issues in
                    % estimating the fraction of time assigned to class k (to
                    % recheck)
                    if ~isempty(PH{i,k})
                        UN(i,k) = probSysState*arvRates(:,i,k)*map_mean(PH{i,k})/S(i);
                    end
            end
        end
    end
    
    for k=1:K
        for i=1:M
            if TN(i,k)>0
                RN(i,k) = QN(i,k)./TN(i,k);
            else
                RN(i,k) = 0;
            end
        end
        CN(k) = NK(k)./XN(k);
    end
    
    QN(isnan(QN))=0;
    CN(isnan(CN))=0;
    RN(isnan(RN))=0;
    UN(isnan(UN))=0;
    XN(isnan(XN))=0;
    TN(isnan(TN))=0;
    
    % now update the routing probabilities in nodes with state-dependent routing
    for k=1:K
        for isf=1:qn.nstateful
            if qn.nodetype(isf) == NodeType.Cache
                TNcache(isf,k) = probSysState*depRates(:,isf,k);
            end
        end
    end
    
    % updates cache actual hit and miss data
    for k=1:K
        for isf=1:qn.nstateful
            if qn.nodetype(isf) == NodeType.Cache
                ind = qn.statefulToNode(isf);
                if length(qnc.varsparam{ind}.hitclass)>=k
                    h = qnc.varsparam{ind}.hitclass(k);
                    m = qnc.varsparam{ind}.missclass(k);
                    qn.varsparam{ind}.actualhitprob(k) = TNcache(isf,h)/sum(TNcache(isf,[h,m]));
                    qn.varsparam{ind}.actualmissprob(k) = TNcache(isf,m)/sum(TNcache(isf,[h,m]));
                end
            end
        end
    end
end

nLabs = length(QN);
QN = cellsum(QN)/nLabs;
UN = cellsum(UN)/nLabs;
RN = cellsum(RN)/nLabs;
TN = cellsum(TN)/nLabs;
CN = cellsum(CN)/nLabs;
XN = cellsum(XN)/nLabs;

for k=1:K
    for isf=1:qnc.nstateful
        if qnc.nodetype(isf) == NodeType.Cache
            ind = qnc.statefulToNode(isf);
            qnc.varsparam{ind}.actualhitprob(k) = 0;
            qnc.varsparam{ind}.actualmissprob(k) = 0;
            for l=1:nLabs
                qntmp = qn{l};                
                if length(qntmp.varsparam{ind}.hitclass)>=k                
                    qnc.varsparam{ind}.actualhitprob(k) = qnc.varsparam{ind}.actualhitprob(k) + (1/nLabs) * qntmp.varsparam{ind}.actualhitprob(k);
                    qnc.varsparam{ind}.actualmissprob(k) = qnc.varsparam{ind}.actualmissprob(k) + (1/nLabs) * qntmp.varsparam{ind}.actualmissprob(k);
                end
            end
        end
    end
end

end
