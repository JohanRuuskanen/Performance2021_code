function estVal = estimateAt(self, node)
switch self.options.method
    case 'ubr' % utilization based regression
        estVal = estimatorUBR(self, node);
    case 'ubo' % utilization based optimization
        estVal = estimatorUBO(self, node);
    case 'erps' % extended regression for PS
        estVal = estimatorERPS(self, node);
    otherwise
        error('Unknown inference method: %s.', self.options.method);
end
estVal = estVal(:)';
% update the model parameters
svcProc = node.getService;
for r=1:self.model.getNumberOfClasses
    svcProc{r}.updateMean(estVal(r));
end
end
