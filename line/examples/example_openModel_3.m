%if ~isoctave(), clearvars -except exampleName testName; end
if ~isoctave(), clearvars -except exampleName; end 
model = JMT2LINE('example_openModel_3.jsimg');

options = Solver.defaultOptions;
options.keep = true;
options.verbose = 1;
options.cutoff = [1,1,0;3,3,0;0,0,3]; % works well with 7
options.seed = 23000;
%options.samples=2e4;

%disp('This example shows the execution of the solver on a 1-class 2-node open model.')
% This part illustrates the execution of different solvers
solver={};
solver{end+1} = SolverCTMC(model,options); % CTMC is infinite on this model
solver{end+1} = SolverJMT(model,options);
solver{end+1} = SolverSSA(model,options);
%solver{end+1} = SolverFluid(model,options);
solver{end+1} = SolverMVA(model,options);
%solver{end+1} = SolverMAM(model,options);
%solver{end+1} = SolverNC(model,options);
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    AvgTable{s} = solver{s}.getAvgTable()
    AvgTable{s}
end
