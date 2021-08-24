function Pnir = getProb(self, node, state)
% PNIR = GETPROB(NODE, STATE)
if ~exist('state','var')
    state = self.model.getState{self.model.getStatefulNodeIndex(node)};
end
T0 = tic;
qn = self.model.getStruct;
% now compute marginal probability
if isa(node,'Node')
    ist = self.model.getStationIndex(node);
else
    ist = node;    
end
qn.state{ist} = state;
if isfield(self.result.Prob,'logNormConstAggr') && isfinite(self.result.Prob.logNormConstAggr)
    [Pnir,lG] = solver_nc_marg(qn, self.options, self.result.Prob.logNormConstAggr);
else
    [Pnir,lG] = solver_nc_marg(qn, self.options);
    self.result.Prob.logNormConstAggr = lG;
end
self.result.('solver') = self.getName();
self.result.Prob.marginal = Pnir;
runtime = toc(T0);
self.result.runtime = runtime;
Pnir = Pnir(ist);
end