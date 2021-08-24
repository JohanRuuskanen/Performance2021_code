function model = gallery_mmap1(MAP)
model = Network('M/MAP/1');
%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.FCFS);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass = OpenClass(model, 'myClass');
source.setArrival(oclass, Exp(1));
queue.setService(oclass, MAP);
%% Block 3: topology
model.link(Network.serialRouting(source,queue,sink));
end