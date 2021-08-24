if ~isoctave(), clearvars -except exampleName; end 
fprintf(1,'This example shows a compact solution of a tandem open queueing network.\n');
D = [10,5; 5,9]; % S(i,r) - mean service time of class r at station i
A = [1,2]/20; % A(r) - arrival rate of class r
Z = [1,2;3,4]; % Z(r)  mean service time of class r at delay station i
solver{1} = SolverMVA(Network.tandemPsInf(A,D,Z));
AvgTable{1} = solver{1}.getAvgTable;
AvgTable{1}

% 1-line version: AvgTable=SolverMVA(Network.tandemPsInf(A,D,Z)).getAvgTable

