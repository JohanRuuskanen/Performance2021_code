function model = gallery_literature_gau12_prob68()
model = Network('GautamBookCRC_Problem68');
%% Block 1: nodes
for i=1:5
    nodes{i} = Queue(model, ['Queue',num2str(i)], SchedStrategy.FCFS);
end
%% Block 2: classes
jobclass{1} = ClosedClass(model, 'Class1', 30, nodes{1});
m = [3,4,5,6,2];
s = [6,2,1,3,2];
for i=1:5
    cx = Coxian.fitMeanAndSCV(m(i),s(i)^2/m(i)^2);    
    nodes{i}.setService(jobclass{1}, cx);
end
%% Block 3: topology
P = model.initRoutingMatrix;
P{1} = [0,0.4,0.6,0,0;
        0,0,0,0.8,0.2;
        0,0,0,0.5,0.5;
        1,0,0,0,0;
        1,0,0,0,0];
model.link(P);    
end
