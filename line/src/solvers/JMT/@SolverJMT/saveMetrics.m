function [simElem, simDoc] = saveMetrics(self, simElem, simDoc)
% [SIMELEM, SIMDOC] = SAVEMETRICS(SIMELEM, SIMDOC)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
numOfPerformanceIndices = length(self.model.perfIndex.Avg);
for j=1:(numOfPerformanceIndices)
    currentPerformanceIndex = self.model.perfIndex.Avg{j,1};
    if currentPerformanceIndex.disabled == 0
        performanceNode = simDoc.createElement('measure');
        performanceNode.setAttribute('alpha', num2str(1 - currentPerformanceIndex.simConfInt,2));
        performanceNode.setAttribute('name', strcat('Performance_', int2str(j)));
        performanceNode.setAttribute('nodeType', 'station');
        performanceNode.setAttribute('precision', num2str(currentPerformanceIndex.simMaxRelErr,2));
        if isempty(currentPerformanceIndex.station)
            performanceNode.setAttribute('referenceNode', '');
        else
            performanceNode.setAttribute('referenceNode', currentPerformanceIndex.station.name);
        end
        performanceNode.setAttribute('referenceUserClass', currentPerformanceIndex.class.name);
        performanceNode.setAttribute('type', currentPerformanceIndex.type);
        performanceNode.setAttribute('verbose', 'false');
        simElem.appendChild(performanceNode);
    end
end
end
