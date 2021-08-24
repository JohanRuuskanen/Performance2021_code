model = Network('model');

source = Source(model,'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model,'Sink');

jobclass = OpenClass(model, 'OpenClass', 0);

source.setArrival(jobclass, Exp(1));
cwd = fileparts(mfilename('fullpath'));
queue.setService(jobclass, Replayer([cwd,filesep,'example_trace.txt']));

model.link(Network.serialRouting(source,queue,sink));

AvgTable = SolverJMT(model,'seed',23000).getAvgTable
