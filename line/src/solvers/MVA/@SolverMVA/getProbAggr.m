function [Pnir,logPnir] = getProbAggr(self, ist)
% [PNIR,LOGPNIR] = GETPROBSTATEAGGR(IST)

if ~exist('ist','var')
    line_error(mfilename,'getProbAggr requires to pass a parameter the station of interest.');
end
if ist > self.model.getNumberOfStations
    line_error(mfilename,'Station number exceeds the number of stations in the model.');
end
if isempty(self.result)
    self.run;
end
Q = self.result.Avg.Q;
qn = self.model.getStruct;
N = qn.njobs;
if all(isfinite(N))
    switch options.method
        case 'exact'
            line_error(mfilename,'Exact marginal state probabilities not available yet in SolverMVA.');
        otherwise
            state = self.model.getState{qn.stationToStateful(ist)};
            [~, nir, ~, ~] = State.toMarginal(qn, ist, state, self.getOptions);
            % Binomial approximation with mean fitted to queue-lengths.
            % Rainer Schmidt, "An approximate MVA ...", PEVA 29:245-254, 1997.
            logPnir = 0;
            for r=1:size(nir,2)
                logPnir = logPnir + nchoosekln(N(r),nir(r));
                logPnir = logPnir + nir(r)*log(Q(ist,r)/N(r));
                logPnir = logPnir + (N(r)-nir(r))*log(1-Q(ist,r)/N(r));
            end
            Pnir = real(exp(logPnir));
    end
else
    line_error(mfilename,'getProbAggr not yet implemented for models with open classes.');
end
end