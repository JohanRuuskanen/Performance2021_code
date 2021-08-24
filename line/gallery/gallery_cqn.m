function model = gallery_cqn(M, useDelay)
model = Network('Single-class CQN');
%% Block 1: nodes
if nargin==1
    useDelay = false;
end
for i=1:M
    node{i} = Queue(model, ['Queue ',num2str(i)], SchedStrategy.PS); 
end
if useDelay
    node{M+1} = DelayStation(model, 'Delay 1');
end
%% Block 2: classes
jobclass{1} = ClosedClass(model, 'Class1', M+3, node{1}, 0);

for i=1:M
    node{i}.setService(jobclass{1}, Exp.fitMean(1.000000+i)); % (Queue 1,Class1)
end
if useDelay
    node{M+1}.setService(jobclass{1}, Exp.fitMean(2.000000)); % (Delay 1,Class1)
end

%% Block 3: topology
P = model.initRoutingMatrix(); % initialize routing matrix 
P{jobclass{1},jobclass{1}} = Network.serialRouting(node);
model.link(P);
end