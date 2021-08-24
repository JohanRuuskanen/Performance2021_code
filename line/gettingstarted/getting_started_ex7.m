model = Network('model');

% Block 1: nodes
node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);

% Block 2: classes
jobclass{1} = ClosedClass(model, 'Class1', 5, node{1}, 0);
node{1}.setService(jobclass{1}, Exp(1.0));
node{2}.setService(jobclass{1}, Exp(0.5));

% Block 3: topology
model.link(Network.serialRouting(node{1},node{2}));

% Block 4: solution
RDfluid = SolverFluid(model).getCdfRespT();
RDsim = SolverJMT(model,'seed',23000,'samples',1e4).getCdfRespT();

%% Plot results
semilogx(RDsim{2,1}(:,2),1-RDsim{2,1}(:,1),'r'); hold all;
semilogx(RDfluid{2,1}(:,2),1-RDfluid{2,1}(:,1),'k--');
legend('jmt-transient','fluid-steady','Location','Best'); 
ylabel('Pr(T > t)'); xlabel('time t');
