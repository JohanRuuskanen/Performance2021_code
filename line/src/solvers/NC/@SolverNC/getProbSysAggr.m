function Pn = getProbSysAggr(self)
% PN = GETPROBSYSSTATEAGGR()

T0 = tic;
qn = self.model.getStruct;
% now compute marginal probability
[Pn,lG] = solver_nc_jointaggr(qn, self.options);
self.result.('solver') = self.getName();
self.result.Prob.logNormConstAggr = lG;
self.result.Prob.joint = Pn;
runtime = toc(T0);
self.result.runtime = runtime;
end