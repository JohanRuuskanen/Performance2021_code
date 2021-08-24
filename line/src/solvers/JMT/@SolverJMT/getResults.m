function [result, parsed] = getResults(self)
% [RESULT, PARSED] = GETRESULTS()

options = self.getOptions;
switch options.method
    case {'jsim','default'}
        [result, parsed] = self.getResultsJSIM;
    otherwise
        [result, parsed] = self.getResultsJMVA;
end
end