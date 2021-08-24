function RD = getCdfRespT(self, R)
% RD = GETCDFRESPT(R)

T0 = tic;
if ~exist('R','var')
    R = self.model.getAvgRespTHandles;
end
qn = self.getStruct;
self.getAvg; % get steady-state solution
options = self.getOptions;
RD = solver_mam_passage_time(qn, qn.ph, options);
runtime = toc(T0);
self.setDistribResults(RD, runtime);
end