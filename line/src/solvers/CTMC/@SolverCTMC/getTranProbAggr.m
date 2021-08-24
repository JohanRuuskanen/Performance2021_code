function [Pi_t, SSnode_a] = getTranProbAggr(self, node)
% [PI_T, SSNODE_A] = GETTRANPROBSTATEAGGR(NODE)

options = self.getOptions;
if isfield(options,'timespan')  && isfinite(options.timespan(2))
    qn = self.getStruct;
    [t,pi_t,~,~,~,~,~,~,~,~,~,SSa] = solver_ctmc_transient_analysis(qn, options);
    jnd = self.model.getNodeIndex(node);
    SSnode_a = SSa(:,(jnd-1)*qn.nclasses+1:jnd*qn.nclasses);
    Pi_t = [t, pi_t];
else
    line_error(mfilename,'getTranProbAggr in SolverCTMC requires to specify a finite timespan T, e.g., SolverCTMC(model,''timespan'',[0,T]).');
end
end