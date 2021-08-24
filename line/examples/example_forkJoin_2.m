model = Network('model');

source = Source(model,'Source');
queue1 = Queue(model,'Queue1',SchedStrategy.PS);
queue2 = Queue(model,'Queue2',SchedStrategy.PS);
fork = Fork(model,'Fork');
fork.setTasksPerLink(2);
join = Join(model,'Join');
sink = Sink(model,'Sink');

jobclass1 = OpenClass(model, 'class1');
jobclass2 = OpenClass(model, 'class2');

source.setArrival(jobclass1, Exp(0.25));
queue1.setService(jobclass1, Exp(1.0));
queue2.setService(jobclass1, Exp(0.75));

source.setArrival(jobclass2, Exp(0.25));
queue1.setService(jobclass2, Immediate());
queue2.setService(jobclass2, Exp(2.0));

P = cellzeros(2,2,6,6);
P{jobclass1,jobclass1}(source,fork) = 1;
P{jobclass1,jobclass1}(fork,queue1) = 1.0;
P{jobclass1,jobclass1}(fork,queue2) = 1.0;
P{jobclass1,jobclass1}(queue1,join) = 1.0;
P{jobclass1,jobclass1}(queue2,join) = 1.0;
P{jobclass1,jobclass1}(join,sink) = 1.0;

P{jobclass2,jobclass2}(source,fork) = 1;
P{jobclass2,jobclass2}(fork,queue1) = 1.0;
P{jobclass2,jobclass2}(fork,queue2) = 1.0;
P{jobclass2,jobclass2}(queue1,join) = 1.0;
P{jobclass2,jobclass2}(queue2,join) = 1.0;
P{jobclass2,jobclass2}(join,sink) = 1.0;

model.link(P);
SolverJMT(model,'keep',true).getAvgNodeTable