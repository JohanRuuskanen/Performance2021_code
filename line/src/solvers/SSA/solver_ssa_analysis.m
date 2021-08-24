function [QN,UN,RN,TN,CN,XN,runtime,tranSysState,tranSync, qnc] = solver_ssa_analysis(qn, options)
% [QN,UN,RN,TN,CN,XN,RUNTIME,TRANSYSSTATE] = SOLVER_SSA_ANALYSIS(QN, OPTIONS)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

M = qn.nstations;    %number of stations
K = qn.nclasses;    %number of classes

S = qn.nservers;
NK = qn.njobs';  % initial population per class
sched = qn.sched;
Tstart = tic;

PH = qn.proc;
tranSysState = [];
tranSync = [];

qnc = qn.copy;
qnc.space = qn.state; % SSA progressively grows this cell array into the simulated state space
switch options.method
    case {'serial.hash','serial.hashed','hashed'}
        options.samples = options.samples + 1;
        [probSysState,SSq,arvRates,depRates] = solver_ssa_hashed(qnc, options);
        qn.space = qnc.space;
        XN = NaN*zeros(1,K);
        UN = NaN*zeros(M,K);
        QN = NaN*zeros(M,K);
        RN = NaN*zeros(M,K);
        TN = NaN*zeros(M,K);
        CN = NaN*zeros(1,K);
        for k=1:K
            refsf = qn.stationToStateful(qn.refstat(k));
            XN(k) = probSysState*arvRates(:,refsf,k);
            for ist=1:M
                isf = qn.stationToStateful(ist);
                TN(ist,k) = probSysState*depRates(:,isf,k);
                QN(ist,k) = probSysState*SSq(:,(ist-1)*K+k);
                switch sched(ist)
                    case SchedStrategy.INF
                        UN(ist,k) = QN(ist,k);
                    otherwise
                        % we use Little's law, otherwise there are issues in
                        % estimating the fraction of time assigned to class k (to
                        % recheck)
                        if ~isempty(PH{ist,k})
                            UN(ist,k) = probSysState*arvRates(:,ist,k)*map_mean(PH{ist,k})/S(ist);
                        end
                end
            end
        end
        
        for k=1:K
            for ist=1:M
                if TN(ist,k)>0
                    RN(ist,k) = QN(ist,k)./TN(ist,k);
                else
                    RN(ist,k) = 0;
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
                    ind = qn.statefulToNode(isf);
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
        qnc = qn;
    case {'default','serial'}
        options.samples = options.samples + 1;
        [probSysState,SSq,arvRates,depRates,tranSysState, tranSync] = solver_ssa(qnc, options);
        qn.space = qnc.space;
        XN = NaN*zeros(1,K);
        UN = NaN*zeros(M,K);
        QN = NaN*zeros(M,K);
        RN = NaN*zeros(M,K);
        TN = NaN*zeros(M,K);
        CN = NaN*zeros(1,K);
        for k=1:K
            refsf = qn.stationToStateful(qn.refstat(k));
            XN(k) = probSysState*arvRates(:,refsf,k);
            for ist=1:M
                isf = qn.stationToStateful(ist);
                TN(ist,k) = probSysState*depRates(:,isf,k);
                QN(ist,k) = probSysState*SSq(:,(ist-1)*K+k);
                switch sched(ist)
                    case SchedStrategy.INF
                        UN(ist,k) = QN(ist,k);
                    otherwise
                        % we use Little's law, otherwise there are issues in
                        % estimating the fraction of time assigned to class k (to
                        % recheck)
                        if ~isempty(PH{ist,k})
                            UN(ist,k) = probSysState*arvRates(:,ist,k)*map_mean(PH{ist,k})/S(ist);
                        end
                end
            end
        end
        
        for k=1:K
            for ist=1:M
                if TN(ist,k)>0
                    RN(ist,k) = QN(ist,k)./TN(ist,k);
                else
                    RN(ist,k) = 0;
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
                    ind = qn.statefulToNode(isf);
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
                        if h==0 | m==0
                            qn.varsparam{ind}.actualhitprob(k) = NaN;
                            qn.varsparam{ind}.actualmissprob(k) = NaN;
                        else
                            qn.varsparam{ind}.actualhitprob(k) = TNcache(isf,h)/sum(TNcache(isf,[h,m]));
                            qn.varsparam{ind}.actualmissprob(k) = TNcache(isf,m)/sum(TNcache(isf,[h,m]));
                        end
                    end
                end
            end
        end
        qnc = qn;
    case {'para','parallel','para.hash','parallel.hash'}
        if isoctave
            line_error(mfilename,'parallel SSA is available only under MATLAB.');
        end
        laboptions = options;
        [XN,UN,QN,RN,TN,CN,qnc] = solver_ssa_analysis_spmd(laboptions, qn, qnc, PH);
end

runtime = toc(Tstart);
end
