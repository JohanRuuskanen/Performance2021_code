function estVal = estimatorUBR(self, node)
% rescale utilization to be mean number of busy serveers
if isfinite(node.getNumberOfServers())
    U = self.getAggrUtil(node);
    if ~isempty(U)
        avgAggrUtil = U.data * node.getNumberOfServers();
    end
end

% obtain per class metrics
isUtilKnown = false(1,self.model.getNumberOfClasses);
for r=1:self.model.getNumberOfClasses
    Ur=self.getUtil(node, self.model.classes{r});
    if isempty(Ur)
        avgUtil{r} = [];
    else
        avgUtil{r} = Ur.data * node.getNumberOfServers();
        isUtilKnown(r) = true;
    end
    avgArvR{r} = self.getArvR(node, self.model.classes{r}).data;
end

% convert cell arrays to scalar
try
    avgA = cell2mat(avgArvR);
    sumUr = sum(cell2mat(avgUtil),2);
    if isempty(sumUr)
        sumUr = 0;
    end
catch me
    switch me.identifier
        case 'MATLAB:catenate:dimensionMismatch'
            error('Sampled metrics have different number of samples, use interpolate() before starting this estimation algorithm.');
    end
end

estVal = zeros(1,self.model.getNumberOfClasses);
% for known pairs arrival rate - per-class utilization
% do separate inference of the demand
for r = find(isUtilKnown)
    estVal(r) = lsqnonneg(avgA(:,r), avgUtil{r});
end

% for the other aggregate classes, do a non-negative
% least squares
estVal(~isUtilKnown) = lsqnonneg(avgA(:,~isUtilKnown), avgAggrUtil - sumUr);

estVal = estVal(:)';
end