function qn = example_randomEnvironment_genqn(rate, N)
%% qn1
qn = Network('qn1');

node{1} = Delay(qn, 'Queue1');
node{2} = Queue(qn, 'Queue2', SchedStrategy.PS);


jobclass{1} = ClosedClass(qn, 'Class1', N, node{1}, 0);

node{1}.setService(jobclass{1}, Exp(rate(1)));
node{2}.setService(jobclass{1}, Exp(rate(2)));

K = 1;
P = cell(K,K);
P{1} = circul(length(node));

qn.link(P);
end