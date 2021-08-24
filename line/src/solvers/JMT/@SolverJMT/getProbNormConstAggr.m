function lNormConst = getProbNormConstAggr(self)
% LNORMCONST = GETPROBNORMCONST()

switch self.options.method
    case {'jmva','jmva.recal','jmva.comom','jmva.ls'}
        self.runAnalysis();
        lNormConst = self.result.Prob.logNormConstAggr;
    otherwise
        lNormConst = NaN; %#ok<NASGU>
        line_error(mfilename,'Selected solver method does not compute normalizing constants. Choose either jmva.recal, jmva.comom, or jmva.ls.');
end
end