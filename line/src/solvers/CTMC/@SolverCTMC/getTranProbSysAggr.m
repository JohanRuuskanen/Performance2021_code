function [Pi_t, SSsysa] = getTranProbSysAggr(self)
% [PI_T, SSSYSA] = GETTRANPROBSYSSTATEAGGR()

options = self.getOptions;
if isfield(options,'timespan')  && isfinite(options.timespan(2))
    qn = self.getStruct;
    [t,pi_t,~,~,~,~,~,~,~,~,~,SSsysa] = solver_ctmc_transient_analysis(qn, options);
    Pi_t = [t, pi_t];
else
    line_error(mfilename,'getTranProbSysAggr in SolverCTMC requires to specify a finite timespan T, e.g., SolverCTMC(model,''timespan'',[0,T]).');
end
end