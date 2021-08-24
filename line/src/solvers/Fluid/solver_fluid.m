function [QN,ymean,QNt,UNt,ymean_t,t,iters,runtime] = solver_fluid(qn, options)
% [QN,YMEAN,QNT,UNT,YMEAN_T,T,ITERS,RUNTIME] = SOLVER_FLUID(QN, OPTIONS)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

M = qn.nstations;    %number of stations
K = qn.nclasses;    %number of classes
N = qn.nclosedjobs;    %population
Mu = qn.mu;
Phi = qn.phi;
PH = qn.proc;
sched = qn.sched;
rt = qn.rt;
S = qn.nservers;
NK = qn.njobs';  %initial population

match = zeros(M,K); % indicates whether a class is served at a station
for i = 1:M
    for k = 1:K
        if isnan(Mu{i,k})
            Mu{i,k} = [];
            Phi{i,k}=[];
        end
        match(i,k) = sum(rt(:, (i-1)*K+k)) > 0;
    end
    %Set number of servers in delay station = population
    if isinf(S(i))
        S(i) = N;
    end
end

%% initialization
max_time = Inf;
Tstart = tic;

%phases = qn.phases;
phases = zeros(M,K);
for i = 1:M
    for k = 1:K
        phases(i,k) = length(Mu{i,k});
    end
end

slowrate = zeros(M,K);
for i = 1:M
    for k = 1:K
        if ~isempty(Mu{i,k})
            slowrate(i,k) = min(Mu{i,k}(:)); %service completion (exit) rates in each phase
        else
            slowrate(i,k) = Inf;
        end
    end
end

%NK(isinf(NK))=1e6;
iters = 0;
ymean = {};
y0 = [];
assigned = zeros(1,K); %number of jobs of each class already assigned
for i = 1:M
    for k = 1:K
        if match(i,k) > 0 % indicates whether a class is served at a station
            if isinf(NK(k))
                if sched(i)==SchedStrategy.EXT
                    toAssign = 1; % open job pool
                else
                    toAssign = 0; % set to zero open jobs everywhere
                end
            else
                toAssign = floor(NK(k)/sum(match(:,k)));
                if sum( match(i+1:end, k) ) == 0 % if this is the last station for this job class
                    toAssign = NK(k) - assigned(k);
                end
            end
            y0 = [y0, toAssign, zeros(1,phases(i,k)-1)]; % this is just for PS
            assigned(k) = assigned(k) + toAssign;
        else
            y0 = [y0, zeros(1,phases(i,k))];
        end
    end
end

if isempty(options.init_sol)
    ymean{iters +1} = y0(:)'; % average state embedded at stage change transitions out of e
    ydefault = y0(:)'; % not used in this case
else
    ymean{iters +1} = options.init_sol(:)';
    ydefault = y0(:)'; % default solution if init_sol fails
end

%% solve ode
[ymean, ymean_t, t, iters] = solver_fluid_iteration(qn, N, Mu, Phi, PH, rt, S, ymean, ydefault, slowrate, Tstart, max_time, options);

runtime = toc(Tstart);
% if options.verbose >= 2
%     if iters==1
%         line_printf('Fluid analysis iteration completed in %0.6f sec [%d iteration]\n',runtime,iters);
%     else
%         line_printf('Fluid analysis iteration completed in %0.6f sec [%d iterations]\n',runtime,iters);
%     end
% end

% this part assumes PS, DPS, GPS scheduling
QN = zeros(M,K);
QNt = cell(M,K);
Qt = cell(M,1);
UNt = cell(M,K);
for i=1:M
    Qt{i} = 0;
    for k = 1:K
        shift = sum(sum(phases(1:i-1,:))) + sum( phases(i,1:k-1) ) + 1;
        QN(i,k) = sum(ymean{end}(shift:shift+phases(i,k)-1));
        QNt{i,k} = sum(ymean_t(:,shift:shift+phases(i,k)-1),2);
        Qt{i} = Qt{i} + QNt{i,k};
        % computes queue length in each node and stage, summing the total
        % number in service and waiting in that station
        % results are weighted with the stat prob of the stage
    end
end

for i=1:M
    if qn.nservers(i) > 0 % not INF
        for k = 1:K
            UNt{i,k} = min(QNt{i,k} / S(i), QNt{i,k} ./ Qt{i}); % if not an infinite server then this is a number between 0 and 1
            UNt{i,k}(isnan(UNt{i,k})) = 0; % fix cases where qlen is 0
        end
    else % infinite server
        for k = 1:K
            UNt{i,k} = QNt{i,k};
        end
    end
end

return
end
