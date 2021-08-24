function [QN,UN,RN,TN,CN,XN,lG,runtime] = solver_mva_cache_analysis(qn, options)
% [Q,U,R,T,C,X,LG,RUNTIME] = SOLVER_MVA_CACHE_ANALYSIS(QN, OPTIONS)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

T0=tic;
QN = []; UN = [];
RN = []; TN = [];
CN = [];
XN = zeros(1,qn.nclasses);
lG = NaN;

source_ist = qn.nodeToStation(qn.nodetype == NodeType.Source);
sourceRate = qn.rates(source_ist,:);
sourceRate(isnan(sourceRate)) = 0;
TN(source_ist,:) = sourceRate;

ch = qn.varsparam{qn.nodetype == NodeType.Cache};

m = ch.cap;
n = ch.nitems;
h = length(m);
u = qn.nclasses;
lambda = zeros(u,n,h);

for v=1:u
    for k=1:n
        for l=1:(h+1)
            if ~isnan(ch.pref{v})
                lambda(v,k,l) = sourceRate(v) * ch.pref{v}(k);
            end
        end
    end
end

R = ch.accost;
gamma = mucache_gamma_lp(lambda,R);

switch options.method
    case 'exact'
        [~,~,pij] = mucache_mva(gamma, m);
        pij = [abs(1-sum(pij,2)),pij];
        missRate = zeros(1,u);
        for v=1:u
            missRate(v) = lambda(v,:,1)*pij(:,1);
        end
    otherwise
        pij = mucache_prob_asy(gamma,m); % FPI method
        missRate = zeros(1,u);
        for v=1:u
            missRate(v) = lambda(v,:,1)*pij(:,1);
        end
end

for r = 1:qn.nclasses
    if length(ch.hitclass)>=r & ch.missclass(r)>0 & ch.hitclass(r)>0
        XN(ch.missclass(r)) = XN(ch.missclass(r)) + missRate(r);
        XN(ch.hitclass(r)) = XN(ch.hitclass(r)) + (sourceRate(r) - missRate(r));
    end
end
runtime=toc(T0);
end
