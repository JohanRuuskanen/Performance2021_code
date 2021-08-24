function [Pi_t, SSnode] = getTranProb(self, node)
% [PI_T, SSNODE] = GETTRANPROBSTATE(NODE)

options = self.getOptions;
if isfield(options,'timespan')  && isfinite(options.timespan(2))
    qn = self.getStruct;
    [t,pi_t,~,~,~,~,~,~,~,~,SS] = solver_ctmc_transient_analysis(qn, options);
    jnd = self.model.getNodeIndex(node);
    shift = 1;
    for isf = 1:qn.nstateful
        len = length(qn.state{isf});
        if qn.statefulToNode(isf) == jnd
            SSnode = SS(:,shift:shift+len-1);
            break;
        end
        shift = shift+len;
    end
    Pi_t = [t, pi_t];
else
    line_error(mfilename,'getTranProb in SolverCTMC requires to specify a finite timespan T, e.g., SolverCTMC(model,''timespan'',[0,T]).');
end
end