function RD = getCdfRespT(self, R)
% RD = GETCDFRESPT(R)

T0 = tic;
if ~exist('R','var')
    R = self.model.getAvgRespTHandles;
end
qn = self.getStruct;
%self.getAvg; % get steady-state solution
options = self.getOptions;
[lambda,D,N,Z,mu,S]= self.model.getProductFormParameters;
fcfsNodes = find(qn.sched(qn.sched ~= SchedStrategy.INF) == SchedStrategy.FCFS);
T = sum(N) * mean(1./qn.rates(fcfsNodes,:));
%tset = [0:T/100000:2.5*T/1000, T/1000:T/1000:T];
tset = logspace(-5,log10(T),100);
rates = qn.rates(qn.sched == SchedStrategy.FCFS,:);
RD = pfqn_stdf(D,N,Z,S,fcfsNodes,rates,tset);
%RD = pfqn_stdf_heur(D,N,Z,S,fcfsNodes,rates,tset);
runtime = toc(T0);
self.setDistribResults(RD, runtime);
end