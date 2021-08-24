if ~isoctave(), clearvars -except exampleName; end 
model = LayeredNetwork('LQN1');

% definition of processors, tasks and entries
P1 = Processor(model, 'P1', 1, SchedStrategy.PS);
T1 = Task(model, 'T1', 1, SchedStrategy.REF).on(P1);
E1 = Entry(model, 'E1').on(T1);

P2 = Processor(model, 'P2', 1, SchedStrategy.PS);
T2 = Task(model, 'T2', 1, SchedStrategy.INF).on(P2);
E2 = Entry(model, 'E2').on(T2);

% definition of activities
T1.setThinkTime(Erlang.fitMeanAndOrder(10,2));

A1 = Activity(model, 'A1', Exp(1.3)).on(T1).boundTo(E1).synchCall(E2,3);
A2 = Activity(model, 'A2', Cox2.fitMeanAndSCV(2.5,10)).on(T2).boundTo(E2).repliesTo(E2);

% instantiate solvers
options = SolverLQNS.defaultOptions;
options.keep = true;
options.verbose = 1;
%options.method = 'lqsim';
%options.samples = 1e4;
lqnssolver = SolverLQNS(model, options);
AvgTableLQNS = lqnssolver.getAvgTable;
AvgTableLQNS

% this method runs the MVA solver in each layer
lnoptions = SolverLN.defaultOptions;
lnoptions.verbose = 1;
lnoptions.seed = 2300;
options = SolverMVA.defaultOptions;
solver{1} = SolverLN(model, @(model) SolverMVA(model, options), lnoptions);
AvgTable{1} = solver{1}.getAvgTable
AvgTable{1}

% this method adapts with the features of each layer
solver{2} = SolverLN(model, @(model) SolverAuto(model, SolverAuto.defaultOptions), lnoptions);
AvgTable{2} = solver{2}.getAvgTable
AvgTable{2}
