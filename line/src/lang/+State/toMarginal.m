function [ni, nir, sir, kir] = toMarginal(qn, ind, state_i, K, Ks, space_buf, space_srv, space_var) %#ok<INUSD>
% [NI, NIR, SIR, KIR] = TOMARGINAL(QN, IND, STATE_I, K, KS, SPACE_BUF, SPACE_SRV, SPACE_VAR) %#OK<INUSD>

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

if ~isa(qn,'NetworkStruct') % the input can be a Network object too
    qn=qn.getStruct();
end

% ind: node index
if ~qn.isstation(ind) && qn.isstateful(ind) % if stateful node
    ni = sum(state_i(1:(end-sum(qn.nvars(ind,:)))));
    nir = state_i(1:(end-sum(qn.nvars(ind,:))));
    sir = nir;
    kir = sir;
    return
end

R = qn.nclasses;
ist = qn.nodeToStation(ind);
%isf = qn.nodeToStateful(ind);

if nargin < 5
    K = qn.phasessz(ist,:);
    Ks = qn.phaseshift(ist,:);
end

isExponential = false;
if max(K)==1
    isExponential = true;
end

if nargin < 8    
    space_var = state_i(:,(end-sum(qn.nvars(ind,:))+1):end); % server stat
    space_srv = state_i(:,(end-sum(K)-sum(qn.nvars(ind,:))+1):(end-sum(qn.nvars(ind,:))));
    space_buf = state_i(:,1:(end-sum(K)-sum(qn.nvars(ind,:))));
end

if isExponential
    sir = space_srv;
    kir = space_srv;
else
    nir = zeros(size(state_i,1),R);
    sir = zeros(size(state_i,1),R); % class-r jobs in service
    kir = zeros(size(state_i,1),R,max(K)); % class-r jobs in service in phase k
    for r=1:R
        for k=1:K(r)
            kir(:,r,k) = space_srv(:,Ks(r)+k);
            sir(:,r) = sir(:,r) + kir(:,r,k);
        end
    end
end
switch qn.schedid(ist)
    case SchedStrategy.ID_INF
        for r=1:R
            nir(:,r) = sir(:,r); % class-r jobs in station
        end
    case SchedStrategy.ID_PS
        for r=1:R
            nir(:,r) = sir(:,r) ; % class-r jobs in station
        end
    case SchedStrategy.ID_EXT
        for r=1:R
            nir(:,r) = Inf;
        end
    case SchedStrategy.ID_FCFS
        for r=1:R
            nir(:,r) = sir(:,r) + sum(space_buf==r,2); % class-r jobs in station
        end
    case SchedStrategy.ID_DPS
        for r=1:R
            nir(:,r) = sir(:,r) ; % class-r jobs in station
        end
    case SchedStrategy.ID_GPS
        for r=1:R
            nir(:,r) = sir(:,r) ; % class-r jobs in station
        end
    case SchedStrategy.ID_HOL
        for r=1:R
            nir(:,r) = sir(:,r) + sum(space_buf==r,2); % class-r jobs in station
        end
    case SchedStrategy.ID_LCFS
        for r=1:R
            nir(:,r) = sir(:,r) + sum(space_buf==r,2); % class-r jobs in station
        end
    case SchedStrategy.ID_SIRO
        for r=1:R
            nir(:,r) = sir(:,r) + space_buf(:,r); % class-r jobs in station
        end
    case SchedStrategy.ID_SEPT
        for r=1:R
            nir(:,r) = sir(:,r) + space_buf(:,r); % class-r jobs in station
        end
    case SchedStrategy.ID_LEPT
        for r=1:R
            nir(:,r) = sir(:,r) + space_buf(:,r); % class-r jobs in station
        end
    otherwise % possibly other stateful nodes
        for r=1:R
            nir(:,r) = sir(:,r) ; % class-r jobs in station
        end
end

for r=1:R
    if isnan(qn.rates(ist,r)) % if disabled
        nir(:,r) = 0;
        for k=1:K(r)
            kir(:,r,k) = 0;
        end
        sir(:,r)=0;
    end
end

ni = sum(nir,2); % total jobs in station
end
