function [QN,UN,RN,TN,CN,XN,InfGen,StateSpace,StateSpaceAggr,EventFiltration,runtime,fname,qnc] = solver_ctmc_analysis(qn, options)
% [QN,UN,RN,TN,CN,XN,INFGEN,STATESPACE,STATESPACEAGGR,EVENTFILTRATION,RUNTIME,FNAME,qn] = SOLVER_CTMC_ANALYSIS(qn, OPTIONS)
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

%if options.remote
%    qn.rtfun = {};
%    qn.lst = {};
%    qn_json = jsonencode(qn);
%    qn = NetworkStruct.fromJSON(qn_json)
    %return
%end

M = qn.nstations;    %number of stations
K = qn.nclasses;    %number of classes
rt = qn.rt;
S = qn.nservers;
NK = qn.njobs';  % initial population per class
sched = qn.sched;

Tstart = tic;
PH = qn.proc;

myP = cell(K,K);
for k = 1:K
    for c = 1:K
        myP{k,c} = zeros(M);
    end
end

for i=1:M
    for j=1:M
        for k = 1:K
            for c = 1:K
                % routing table for each class
                myP{k,c}(i,j) = rt((i-1)*K+k,(j-1)*K+c);
            end
        end
    end
end

if any(qn.nodetype == NodeType.Cache)
    options.hide_immediate = false;
end

[InfGen,StateSpace,StateSpaceAggr,EventFiltration,arvRates,depRates,qn] = solver_ctmc(qn, options); % qn is updated with the state space

% if the initial state does not reflect the final state of the state
% vectors, attempt to correct it
for isf=1:qn.nstateful
    if size(qn.state{isf},2) < size(qn.space{isf},2)
        row = matchrow(qn.space{isf}(:,end-length(qn.state{isf})+1:end),qn.state{isf});
        if row > 0
            qn.state{isf} = qn.space{isf}(row,:);
        end
    end
end
qnc = qn;

if options.keep
    fname = tempname;
    save([fname,'.mat'],'InfGen','StateSpace','StateSpaceAggr','EventFiltration')
    line_printf('\nCTMC infinitesimal generator and state space saved in: ');
    line_printf([fname, '.mat'])
else
    fname = '';
end

wset = 1:length(InfGen);

% for caches we keep the immediate transitions to give hit/miss rates
[pi, ~, nConnComp, connComp] = ctmc_solve(InfGen, options);

if any(isnan(pi))
    if nConnComp > 1
        % the matrix was reducible
        initState = matchrow(StateSpace, cell2mat(qn.state'));
        % determine the weakly connected component associated to the initial state
        wset = find(connComp == connComp(initState));
        pi = ctmc_solve(InfGen(wset, wset), options);
        InfGen = InfGen(wset, wset);
        StateSpace = StateSpace(wset,:);
    end
end

pi(pi<1e-14)=0;
pi = pi/sum(pi);

XN = NaN*zeros(1,K);
UN = NaN*zeros(M,K);
QN = NaN*zeros(M,K);
RN = NaN*zeros(M,K);
TN = NaN*zeros(M,K);
CN = NaN*zeros(1,K);

for k=1:K
    refsf = qn.stationToStateful(qn.refstat(k));
    XN(k) = pi*arvRates(wset,refsf,k);
    for i=1:M
        isf = qn.stationToStateful(i);
        TN(i,k) = pi*depRates(wset,isf,k);
        QN(i,k) = pi*StateSpaceAggr(wset,(i-1)*K+k);
        switch sched(i)
            case SchedStrategy.INF
                UN(i,k) = QN(i,k);
            otherwise
                if ~isempty(PH{i,k})
                    UN(i,k) = pi*arvRates(wset,isf,k)*map_mean(PH{i,k})/S(i);
                end
        end
    end
end

for k=1:K
    for i=1:M
        if TN(i,k)>0
            RN(i,k) = QN(i,k)./TN(i,k);
        else
            RN(i,k)=0;
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

runtime = toc(Tstart);

% now update the routing probabilities in nodes with state-dependent routing
for k=1:K
    for isf=1:qn.nstateful
        if qn.nodetype(isf) == NodeType.Cache            
            TNcache(isf,k) = pi*depRates(wset,isf,k);
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
