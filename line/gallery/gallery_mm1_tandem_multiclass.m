function model = gallery_mm1_tandem_multiclass
model = Network('M[2]/M[2]/1 -> -/M[2]/1');
%% Block 1: nodes
source = Source(model, 'Source');
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS);
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass1 = OpenClass(model, 'myClass1');
source.setArrival(oclass1, Exp(1));
queue1.setService(oclass1, Exp(4));
queue2.setService(oclass1, Exp(6));

oclass2 = OpenClass(model, 'myClass2');
source.setArrival(oclass2, Exp(0.5));
queue1.setService(oclass2, Exp(2));
queue2.setService(oclass2, Exp(6));
%% Block 3: topology
P = model.initRoutingMatrix;
P{1} = Network.serialRouting(source,queue1,queue2,sink);
P{2} = Network.serialRouting(source,queue1,queue2,sink);
model.link(P);
end