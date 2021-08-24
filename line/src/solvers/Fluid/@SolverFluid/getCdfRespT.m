function RD = getCdfRespT(self, R)
% RD = GETCDFRESPT(R)

T0 = tic;
if ~exist('R','var')
    R = self.model.getAvgRespTHandles;
end
qn = self.getStruct;
self.getAvg; % get steady-state solution
options = self.getOptions;
options.init_sol = self.result.solverSpecific.odeStateVec;
RD = solver_fluid_passage_time(qn, options);
runtime = toc(T0);
self.setDistribResults(RD, runtime);
end