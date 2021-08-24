if ~isoctave(), clearvars -except exampleName; end
N = 1;
M = 2;

E = 2;
envModel = Env('MyEnv');
envName = {'Stage1', 'Stage2'};
envType = {'UP', 'DOWN'};

rate = zeros(M,E); rate(M,1:E)=(1:E); rate(1,1:E)=(E:-1:1);
envSubModel = {example_randomEnvironment_genqn(rate(:,1),N), example_randomEnvironment_genqn(rate(:,2),N)};
for e=1:E
    envModel.addStage(envName{e}, envType{e}, envSubModel{e});
end

envRates = [0,1; 0.5,0.5];
for e=1:E
    for h=1:E
        if envRates(e,h)>0
            envModel.addTransition(envName{e}, envName{h}, Exp(envRates(e,h)));
        end
    end
end

%
fprintf(1,'The metasolver considers an environment with 2 stages and a queueing network with 2 stations.\n')
fprintf(1,'Every time the stage changes, the queueing network will modify the service rates of the stations.\n')

%options.iter_tol = 1e-5;
options = Solver.defaultOptions;
options.timespan = [0,Inf];
options.iter_max = 100;
options.iter_tol = 0.01;
options.method = 'default';
options.verbose = true;

envModel.getStageTable

sfoptions = SolverFluid.defaultOptions;
sfoptions.timespan = [0,1e3];
sfoptions.verbose = false;
envSolver = SolverEnv(envModel,@(model) SolverFluid(model, sfoptions),options);
[QN,UN,TN] = envSolver.getAvg();
AvgTable = envSolver.getAvgTable()

% scoptions = SolverJMT.defaultOptions;
% scoptions.timespan = [0,1e3];
% scoptions.verbose = false;
% envSolver = SolverEnv(envModel,@(model) SolverJMT(model, scoptions),options);
% [QNc,UNc,TNc] = envSolver.getAvg();
% AvgTableC = envSolver.getAvgTable()

