function [ode_h,q_indices] = solver_fluid_odes(N, Mu, Phi, PH, P, nservers, sched, schedparam, options)
% [ODE_H,Q_INDICES] = SOLVER_FLUID_ODES(N, MU, PHI, PH, P, NSERVERS, SCHED, SCHEDPARAM)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

M = length(nservers);    % number of stations
K = size(Mu,2);   % number of classes
w = ones(M,K);

enabled = false(M,K); % indicates whether a class is served at a station
for i = 1:M
    for c = 1:K
        %enabled(i,c) = sum( P(:,(i-1)*K+c) ) > 0;
        %changed to consider rows instead of colums, for response time
        %analysis (first articial class never returns to delay node)
        enabled(i,c) = sum( P((i-1)*K+c,:) ) > 0;
    end
end

q_indices = zeros(M,K);
Kic = zeros(M,K);
cumsum = 1;
for i = 1 : M
    for c = 1:K
        q_indices(i,c) = cumsum;
        if isnan(Mu{i,c})
            numphases = 0;
        else
            numphases = length(Mu{i,c});
        end
        Kic(i,c) = numphases;
        cumsum = cumsum + numphases;
    end
end
% to speed up convert sched strings in numerical values
sched_id = zeros(1,M);
for i = 1 : M
    sched_id(i) = SchedStrategy.toId(sched(i));
    switch sched(i) % source
        case SchedStrategy.DPS
            w(i,:) = schedparam(i,:);
    end
end

%[rateBase, eventIdx] = ode_hybrid_rate_base(Phi, Mu, PH, M, K, enabled, q_indices, P, Kic, sched_id, all_jumps);

%% define ODE system to be returned
switch options.method
    case {'default','stateindep'}
        % determine all the jumps, and saves them for later use
        all_jumps = ode_jumps_new(M, K, enabled, q_indices, P, Kic);
        % determines a vector with the fixed part of the rates,
        % and defines the indexes that correspond to the events that occur
        [rateBase, eventIdx] = ode_rate_base(Phi, Mu, PH, M, K, enabled, q_indices, P, Kic, sched_id, all_jumps);
        ode_si_h = @(t,x) all_jumps * ode_rates_stateindep(x, M, K, q_indices, Kic, nservers, w, sched_id, rateBase, eventIdx);
        ode_h = ode_si_h;
    case 'statedep'
        ode_sd_h = @(t,x) ode_statedep(x, Phi, Mu, PH, M, K, enabled, q_indices, P, Kic, nservers, w, sched_id);
        ode_h = ode_sd_h;
        %ode_sd_h = @(t,x) ode_hybrid(x, Phi, Mu, PH, M, K, enabled, q_indices, P, Kic, nservers, w, sched_id, all_jumps, rateBase, eventIdx);
end
end