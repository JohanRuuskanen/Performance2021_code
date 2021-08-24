function [ymean, ymean_t, t, iter] = solver_fluid_iteration(qn, N, Mu, Phi, PH, P, S, ymean, ydefault, slowrate, Tstart, max_time, options)
% [YMEAN, YMEAN_T, T, ITER] = SOLVER_FLUID_ITERATION(QN, N, MU, PHI, PH, P, S, YMEAN, YDEFAULT, SLOWRATE, TSTART, MAX_TIME, OPTIONS)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

iter_max = options.iter_max;
verbose = options.verbose;
tol = options.tol;
iter_tol = options.iter_tol;
stiff = options.stiff;
timespan = options.timespan;

goon = true; % max stiff solver
iter = 0;
allt=[];
ally =[];
lastmsg = '';
t=[];
ymean_t=[];
% heuristic to select stiff or non-stiff ODE solver
nonZeroRates = slowrate(:);
nonZeroRates = nonZeroRates(nonZeroRates>0);

nonZeroRates = nonZeroRates(isfinite(nonZeroRates));
rategap = log10(max(nonZeroRates)/min(nonZeroRates)); % if the max rate is Distrib.InfRate and the min is 1, then rategap = 6

% init ode
[ode_h, ~] = solver_fluid_odes(N, Mu, Phi, PH, P, S, qn.sched, qn.schedparam, options);

T0 = timespan(1);
%opt = odeset();
%opt = odeset('AbsTol', min(10^(-rategap),1e-4), 'RelTol', 1e-3, 'NonNegative', 1:length(y0));
odeopt = odeset('AbsTol', tol, 'RelTol', tol, 'NonNegative', 1:length(ymean{1}));
T = 0;

while (isfinite(timespan(2)) && T < timespan(2)) || (goon && iter < iter_max)
    iter = iter + 1;
    if toc(Tstart) > max_time
        goon = false;
        break;
    end
    
    % determine entry state vector in e
    y0 = ymean{iter-1 +1};
    
    if iter == 1 % first iteration
        T = min(timespan(2),abs(10/min(nonZeroRates))); % solve ode until T = 1 event with slowest exit rate
    else
        T = min(timespan(2),abs(10*iter/min(nonZeroRates)));
    end
    trange = [T0, T];
    
    try
        if stiff
            [t_iter, ymean_t_iter] = solveodestiff(ode_h, trange, y0, odeopt, options);
        else
            [t_iter, ymean_t_iter] = solveode(ode_h, trange, y0, odeopt, options);
        end
    catch me
        line_printf('\nThe initial point is invalid, Fluid solver switching to default initialization.');
        odeopt = odeset('AbsTol', tol, 'RelTol', tol, 'NonNegative', 1:length(ydefault));
        [t_iter, ymean_t_iter] = solveode(ode_h, trange, ydefault, odeopt, options);
    end
    
    ymean_t(end+1:end+size(ymean_t_iter,1),:) = ymean_t_iter;
    t(end+1:end+size(t_iter,1),:) = t_iter;
    ymean{iter +1} = ymean_t(end,:);
    movedMassRatio = norm(ymean{iter +1} - ymean{iter-1 +1}, 1) / 2 / sum(ymean{iter-1 +1});
    T0  = T; % for next iteration
    
    if nargout>3
        if isempty(allt)
            allt = t;
        else
            allt = [allt; allt(end)+t];
        end
        ally = [ally; ymean_t];
    end
    
    % check termination condition
    
    if verbose > 0
        llmsg = length(lastmsg);
        if llmsg>0
            for ib=1:llmsg
                line_printf('\b');
            end
        end
    end
    
    if movedMassRatio < iter_tol  % converged
        % stop only if this is not a transient analysis, in which case keep
        % going until the specified end time
    end
    
    if T >= timespan(2)
        goon = false;
    end
end

end
