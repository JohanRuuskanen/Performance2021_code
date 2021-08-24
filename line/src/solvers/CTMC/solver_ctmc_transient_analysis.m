function [t,pit,QNt,UNt,RNt,TNt,CNt,XNt,InfGen,StateSpace,StateSpaceAggr,EventFiltration,runtime,fname] = solver_ctmc_transient_analysis(qn, options)
% [T,PIT,QNT,UNT,RNT,TNT,CNT,XNT,INFGEN,STATESPACE,STATESPACEAGGR,EVENTFILTRATION,RUNTIME,FNAME] = SOLVER_CTMC_TRANSIENT_ANALYSIS(QN, OPTIONS)
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

RNt=[]; CNt=[];  XNt=[];

M = qn.nstations;    %number of stations
K = qn.nclasses;    %number of classes
fname = '';
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

[InfGen,StateSpace,StateSpaceAggr,EventFiltration,~,depRates,qn] = solver_ctmc(qn, options); % qn is updated with the state space


if options.keep
    fname = tempname;
    save([fname,'.mat'],'InfGen','StateSpace','StateSpaceAggr','EventFiltration')
    line_printf('\nCTMC infinitesimal generator and state space saved in: ');
    line_printf([fname, '.mat'])
end

state = [];
for i=1:qn.nnodes
    if qn.isstateful(i)
        isf = qn.nodeToStateful(i);
        state = [state,zeros(1,size(qn.space{isf},2)-length(qn.state{isf})),qn.state{isf}];
    end
end
pi0 = zeros(1,length(InfGen));

state0 = matchrow(StateSpace, state);
if state0 == -1
    line_error(mfilename,'Initial state not contained in the state space.');   
%     state0 = matchrow(StateSpace, round(state));
%     state = round(state);
%     if state0 == -1
%         line_error(mfilename,'Cannot recover - CTMC stopping');
%     end
end
pi0(state0) = 1; % find initial state and set it to probability 1

%if options.timespan(1) == options.timespan(2)
%    pit = ctmc_uniformization(pi0,Q,options.timespan(1));
%    t = options.timespan(1);
%else
[pit,t] = ctmc_transient(InfGen,pi0,options.timespan(1),options.timespan(2),options.stiff);
%end
pit(pit<1e-14)=0;

QNt = cell(M,K);
UNt = cell(M,K);
%XNt = cell(1,K);
TNt = cell(M,K);

if t(1) == 0
    t(1) = 1e-8;
end
for k=1:K
    %    XNt(k) = pi*arvRates(:,qn.refstat(k),k);
    for i=1:M        
        %occupancy_t = cumsum(pit.*[0;diff(t)],1)./t;        
        occupancy_t = pit;
        TNt{i,k} = occupancy_t*depRates(:,i,k);
        qlenAt_t = pit*StateSpaceAggr(:,(i-1)*K+k);
        %QNt{i,k} = cumsum(qlenAt_t.*[0;diff(t)])./t;
        QNt{i,k} = qlenAt_t;
        switch sched(i)
            case SchedStrategy.INF
                UNt{i,k} = QNt{i,k};
            case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.SIRO, SchedStrategy.SEPT, SchedStrategy.LEPT, SchedStrategy.SJF}
                if ~isempty(PH{i,k})
                    UNt{i,k} = occupancy_t*min(StateSpaceAggr(:,(i-1)*K+k),S(i))/S(i);
                end
            case SchedStrategy.PS
                uik = min(StateSpaceAggr(:,(i-1)*K+k),S(i)) .* StateSpaceAggr(:,(i-1)*K+k) ./ sum(StateSpaceAggr(:,((i-1)*K+1):(i*K)),2);
                uik(isnan(uik))=0;
                utilAt_t = pit * uik / S(i);
                %UNt{i,k} = cumsum(utilAt_t.*[0;diff(t)])./t;
                UNt{i,k} = utilAt_t;
            case SchedStrategy.DPS
                w = qn.schedparam(i,:);
                nik = S(i) * w(k) * StateSpaceAggr(:,(i-1)*K+k) ./ sum(repmat(w,size(StateSpaceAggr,1),1).*StateSpaceAggr(:,((i-1)*K+1):(i*K)),2);
                nik(isnan(nik))=0;
                UNt{i,k} = occupancy_t*nik;
            otherwise
                if ~isempty(PH{i,k})
                    ind = qn.stationToNode(i);
                    line_warning(mfilename,'Transient utilization not support yet for station %s, returning an approximation.',qn.nodenames{ind});
                    UNt{i,k} = occupancy_t*min(StateSpaceAggr(:,(i-1)*K+k),S(i))/S(i);
                end
        end
    end
end
runtime = toc(Tstart);

if options.verbose > 0
    line_printf('\nCTMC analysis completed. Runtime: %f seconds.\n',runtime);
end
end
