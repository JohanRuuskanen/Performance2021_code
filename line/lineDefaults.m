function options = lineDefaults(solverName)
if nargin < 1
    solverName = 'Solver'; % global options unless overridden by a solver
end
%% Solver default options
options = struct();
options.cache = true;
options.cutoff = Inf;
options.config = {};
options.force = false;
options.hide_immediate = true; % hide immediate transitions if possible
options.init_sol = [];
options.iter_max = 10;
options.iter_tol = 1e-4; % convergence tolerance to stop iterations
options.tol = 1e-4; % tolerance for all other uses
options.keep = false;
options.method = 'default';
options.remote = false;
options.remote_endpoint = '127.0.0.1';

odesfun = struct();
odesfun.fastOdeSolver = Solver.fastOdeSolver;
odesfun.accurateOdeSolver = Solver.accurateOdeSolver;
odesfun.fastStiffOdeSolver = Solver.fastStiffOdeSolver;
odesfun.accurateStiffOdeSolver = Solver.accurateStiffOdeSolver;
options.odesolvers = odesfun;

if isoctave
    options.samples = 5e3;
else
    options.samples = 1e4;
end
options.seed = randi([1,1e6]);
options.stiff = true;
options.timespan = [Inf,Inf];
options.verbose = 1;

%% Solver-specific defaults
switch solverName
    case 'CTMC'
        options.timespan = [Inf,Inf];
    case 'Ensemble' % Env
        options.method = 'default';
        options.init_sol = [];
        options.iter_max = 100;
        options.iter_tol = 1e-4;
        options.tol = 1e-4;
        options.verbose = 0;    
    case 'Fluid'
        options = Solver.defaultOptions();
        options.iter_max = 50;
        options.stiff = true;
        options.timespan = [0,Inf];
    case 'JMT'
        % use default
    case 'LN'
        options = EnsembleSolver.defaultOptions();
        options.timespan = [Inf,Inf];
        options.keep = false;
        options.verbose = true;
        options.iter_tol = 5e-2;
        options.iter_max = 100;
    case 'LQNS'
        options = EnsembleSolver.defaultOptions();
        options.timespan = [Inf,Inf];
        options.keep = false;
    case 'MAM'
        options.timespan = [Inf,Inf];
    case 'NC'        
        options.samples = 1e5;
        options.timespan = [Inf,Inf];
    case 'SSA'
        options.timespan = [0,Inf];
        options.verbose = true;
end
end