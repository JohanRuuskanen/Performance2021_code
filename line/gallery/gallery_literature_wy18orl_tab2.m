function model = gallery_literature_wy18orl_tab2()
% Model from "The Advantage of Indices of Dispersion in Queueing
% Approximations", by Ward Whitt and Wei You, 2018, Table 2.
% Expected waiting time (excluding job in service)
% Sim= [3.3600    2.3200    1.9600    1.7700    1.6400    1.5600    1.4900 1.4400   29.2000];
model = Network('WhittYou2018ORL_Table2');
%% Block 1: nodes
node{1} = Source(model, 'mySource');
for i=1:9
    node{1+i} = Queue(model, ['Queue',num2str(i)], SchedStrategy.FCFS);
end
node{end+1} = Sink(model, 'mySink');
%% Block 2: classes
oclass = OpenClass(model, 'myClass');
node{1}.setArrival(oclass, HyperExp.fitMeanAndSCVBalanced(1, 8));
means = [0.6*ones(1,8), 0.9]; 
for i=1:9
    node{1+i}.setService(oclass, Exp.fitMean(means(i)));
end
%% Block 3: topology
model.link(Network.serialRouting(node));
end
