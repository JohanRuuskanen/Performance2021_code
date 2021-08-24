function Pnir = getProbAggr(self, node, state_a)
% PNIR = GETPROBAGGR(NODE, STATE_A)

if ~exist('state_a','var')
    state_a = self.model.getState{self.model.getStatefulNodeIndex(node)};
end
T0 = tic;
qn = self.model.getStruct;
% now compute marginal probability
ist = self.model.getStationIndex(node);
qn.state{ist} = state_a;

self.result.('solver') = self.getName();
if isfield(self.result,'Prob') && isfield(self.result.Prob,'logNormConstAggr') && isfinite(self.result.Prob.logNormConstAggr)
    [Pnir,lG] = solver_nc_margaggr(qn, self.options, self.result.Prob.logNormConstAggr);
else
    [Pnir,lG] = solver_nc_margaggr(qn, self.options);
    self.result.Prob.logNormConstAggr = lG;
end
self.result.Prob.marginal = Pnir;
runtime = toc(T0);
self.result.runtime = runtime;
Pnir = Pnir(ist);
end