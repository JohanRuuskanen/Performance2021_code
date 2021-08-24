function [logNormConst] = getProbNormConstAggr(self)
% [LOGNORMCONST] = GETPROBNORMCONST()

self.runAnalysis();
logNormConst = self.result.Prob.logNormConstAggr;
end