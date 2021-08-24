if ~isoctave(), clearvars -except exampleName testName; end
% this is a difficult sparse LDMX model
model = Network('model');

node{1} = Queue(model, 'Queue1', SchedStrategy.PS);
node{1}.setNumberOfServers(1);
node{2} = Queue(model, 'Queue2', SchedStrategy.PS);
node{2}.setNumberOfServers(2);
node{3} = Queue(model, 'Queue3', SchedStrategy.PS);
node{3}.setNumberOfServers(3);
node{4} = Queue(model, 'Queue4', SchedStrategy.PS); % only closed classes
node{4}.setNumberOfServers(4);
node{5} = Queue(model, 'Queue5', SchedStrategy.PS); % only open classes
node{5}.setNumberOfServers(5);
source = Source(model,'Source');
sink = Sink(model,'Sink');

jobclass{1} = ClosedClass(model, 'Class1', 3, node{1}, 0);
jobclass{2} = OpenClass(model, 'Class2', 0);

for i=1:5
    node{i}.setService(jobclass{1}, Exp(i));
    node{i}.setService(jobclass{2}, Exp(sqrt(i)));
end

source.setArrival(jobclass{2}, Exp(0.3));

M = model.getNumberOfStations();
K = model.getNumberOfClasses();

P = model.initRoutingMatrix;
P{1,1} = zeros(M+1);
P{1,1} = Network.serialRouting(node{1},node{2},node{3},node{4});
P{1,2} = zeros(M+1);
P{2,1} = zeros(M+1);
P{2,2} = Network.serialRouting(source,node{1},node{2},node{3},node{5},sink);

model.link(P);
%%
options = Solver.defaultOptions;
options.keep=false;
options.verbose=1;
options.cutoff = 3;
options.seed = 23000;
%options.samples=2e4;
optionssa = options; optionssa.cutoff = Inf;

disp('This example shows the execution of the solver on a 2-class mixed model with 5 multi-server nodes.')
% This part illustrates the execution of different solvers
solver={};
solver{end+1} = SolverCTMC(model,options); % CTMC is infinite on this model
solver{end+1} = SolverJMT(model,'samples',1e5);
solver{end+1} = SolverSSA(model,options);
%solver{end+1} = SolverFluid(model,options);
solver{end+1} = SolverMVA(model,'exact');
%solver{end+1} = SolverNC(model,options);
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}    
end
