function RTret = solver_fluid_passage_time(qn, options)
% RTRET = SOLVER_FLUID_PASSAGE_TIME(QN, OPTIONS)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

iter_max = options.iter_max;
tol = options.tol;
y0 = options.init_sol;
stiff = options.stiff;
T0 = 0;
M = qn.nstations;    %number of stations
K = qn.nclasses;    %number of classes
N = qn.nclosedjobs;    %population
Lambda = qn.mu;
Pi = qn.phi;
PH = qn.proc;
rt = qn.rt;
S = qn.nservers;

for j = 1:M
    %Set number of servers in delay station = population
    if isinf(S(j))
        S(j) = N;
    end
end

%% initialization
phases = qn.phases;
slowrate = zeros(M,K);
for j = 1:M
    for k = 1:K
        slowrate(j,k) = Inf;
        slowrate(j,k) = min(slowrate(j,k),min(Lambda{j,k}(:))); %service completion (exit) rates in each phase
    end
end

%% response time analysis - starting from fixed point found
%stiff = 1;
chains = qn.chains;
nChains = size(chains,1);

RT = [];
for i = 1:qn.nstations
    for k = 1:nChains %once for each chain
        idxClassesInChain = find(chains(k,:)==1);
        for c = idxClassesInChain
            if phases(i,c) > 0
                [Kc, ode_h_c, y0_c, phases_c, fluid_c] = generate_odes_passagetime(k,c);
                
                % determine max integration time
                nonZeroRates = slowrate(:);
                nonZeroRates = nonZeroRates( nonZeroRates >0 );
                T = abs(100/min(nonZeroRates)); % solve ode until T = 100 events with slowest exit rate
                
                % indices of new classes at station i
                idxN = [];
                for j = i % this used to be at all stations
                    idxN = [idxN sum(sum(phases_c(1:j-1,: ) )) + sum(phases_c(j,1:K)) + [1:sum(phases_c(j,K+1:Kc))] ];
                end
                
                %% ODE analysis
                fullt = [];
                fully = [];
                iter = 1;
                finished = 0;
                tref = 0;
                odeopt = odeset('AbsTol', tol, 'RelTol', tol, 'NonNegative', 1:length(y0_c));
                while iter <= iter_max && finished == 0
                    trange = [T0, T];
                    % solve ode - ymean_t_iter is the transient solution in stage e
                    try
                        if stiff
                            [t_iter, ymean_t_iter] = solveodestiff(ode_h_c, trange, y0_c, odeopt, options);
                        else
                            [t_iter, ymean_t_iter] = solveode(ode_h_c, trange, y0_c, odeopt, options);
                        end
                    catch me
                        line_printf('\nODE solver failed. Fluid solver switching to default initialization.');
                        odeopt = odeset('AbsTol', tol, 'RelTol', tol, 'NonNegative', 1:length(y0_c));
                        try
                            [t_iter, ymean_t_iter] = solveode(ode_h_c, trange, y0_c, odeopt, options);
                        catch
                            keyboard
                        end
                    end
                    %%%
                    iter = iter + 1;
                    fullt = [fullt; t_iter+tref];
                    fully = [fully; ymean_t_iter];
                    if sum(ymean_t_iter(end,idxN )) < 10e-10
                        finished = 1;
                    end
                    tref = tref + t_iter(end);
                    y0_c = ymean_t_iter(end,:);
                end
                
                % retrieve response time CDF for class k
                RT{i,c,1} = fullt;
                if fluid_c > 0
                    RT{i,c,2} = 1 - sum(fully(:,idxN ),2)/fluid_c;
                else
                    RT{i,c,2} = ones(size(fullt));
                end
                if iter > iter_max
                    line_warning(mfilename,'Maximum number of iterations reached when computing the response time distribution. Response time distributions may be inaccurate. Increase option.iter_max (currently at %s).',num2str(iter_max));
                end
            end
        end
    end
end

RTret = {};
for i=1:qn.nstations
    for c=1:qn.nclasses
        RTret{i,c} = [RT{i,c,2},RT{i,c,1}];
    end
end
return

    function [Kc, ode_h_c, y0_c, phases_c, fluid_c] = generate_odes_passagetime(k,c)
        Kc = K + 1;  % add a single new class
        numTranClasses = Kc - K;
        idxTranCl = zeros(1,K); % indices of the transient class corresponding to each class in the original model for chain k
        idxTranCl(chains(k,:)==1) =  K+1:Kc; % this is just K+1 since we are adding a single new class, but this format may be generalizable
        newLambda = cell(M,Kc);
        new_pi = cell(M,Kc);
        new_rt = zeros(M*Kc, M*Kc); % new routing table
        new_proc = PH;
        
        for j=1:qn.nstations
            % service rates
            newLambda(j,1:K) = Lambda(j,:);
            newLambda(j,K+1) = Lambda(j,c);
            
            % completion probabilities
            new_pi(j,1:K) = Pi(j,:);
            new_pi(j,K+1) = Pi(j,c);
            
            % phd distribution
            for r=1:nChains
                new_proc{j,r} = PH{j,r};
            end
            new_proc{j,K+1} = PH{j,c};
        end
        
        % routing/switching probabilities
        % among basic classes
        for l = 1:K
            for m = 1:K
                new_rt(l:Kc:end,m:Kc:end) = rt(l:K:end,m:K:end);
            end
        end
        
        % copy routing table from the original to the transient classes (forward)
        for l = 1:numTranClasses
            for m = 1:numTranClasses
                if sum(sum(rt(c:K:end,idxClassesInChain(m):K:end))) > 0
                    new_rt(K+l:Kc:end,K+m:Kc:end) = rt(c:K:end,idxClassesInChain(m):K:end);
                end
            end
        end
        
        %phases of transient classes
        phases_c = zeros(M,Kc);
        phases_c(:,1:K) = phases;
        phases_c(:,K+1) = phases(:,c);
        
        % identify classes in chain that complete
        completingClassesInChain = c;
        
        %determine final classes (leaves in the class graph)
        for s = completingClassesInChain' % for each completing class
            %routing matrix from a transient class that completes is diverted back into the original classes
            for l = idxClassesInChain
                for j = 1:qn.nstations
                    % return fluid to original class
                    new_rt((i-1)*Kc+idxTranCl(c), (j-1)*Kc+l) = rt((i-1)*K+c, (j-1)*K+l);
                    % delete corresponding transition among transient classes
                    new_rt((i-1)*Kc+idxTranCl(c), (j-1)*Kc+idxTranCl(l)) = 0;
                end
            end
        end
        
        % setup the ODEs for the new QN
        %        options.method  = 'statedep'; % default doesn't seem to work in some models
        [ode_h_c, ~] = solver_fluid_odes(N, reshape({newLambda{:,:}},M,Kc), reshape({new_pi{:,:}},M,Kc), new_proc, new_rt, S, qn.sched, qn.schedparam, options);
        
        % setup initial point
        y0_c = zeros(1, sum(sum(phases_c(:,:))));
        fluid_c = 0;
        for j = 1:qn.nstations
            for l = 1:qn.nclasses
                idxNew_jl = sum(sum(phases_c(1:j-1,:))) + sum(phases_c(j,1:l-1));
                idxNew_jt = sum(sum(phases_c(1:j-1,:))) + sum(phases_c(j,1:idxTranCl(l)-1));
                idx_jl = sum(sum(phases(1:j-1,:))) + sum(phases(j,1:l-1));
                if i == j && l==c
                    y0_c( idxNew_jt + 1 ) = sum(y0(idx_jl+1:idx_jl + phases(j,l))); % mass in phases all moved back into phase 1
                    fluid_c = fluid_c + sum(y0(idx_jl+1:idx_jl + phases(j,l)));
                else % leave mass as it is
                    y0_c( idxNew_jl + 1: idxNew_jl + phases_c(j,l)  ) = y0(idx_jl+1:idx_jl + phases(j,l));
                end
            end
        end
    end

end