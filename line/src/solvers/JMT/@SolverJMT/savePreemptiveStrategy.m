function [simDoc, section] = savePreemptiveStrategy(self, simDoc, section, currentNode)
% [SIMDOC, SECTION] = SAVEPREEMPTIVESTRATEGY(SIMDOC, SECTION, CURRENTNODE)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
visitsNode = simDoc.createElement('parameter');
visitsNode.setAttribute('array', 'true');
visitsNode.setAttribute('classPath', 'jmt.engine.NetStrategies.PSStrategy');
visitsNode.setAttribute('name', 'PSStrategy');

numOfClasses = length(self.model.classes);
for i=1:(numOfClasses)
    currentClass = self.model.classes{i,1};
    
    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(currentClass.name));
    visitsNode.appendChild(refClassNode);
    
    subParameterNode = simDoc.createElement('subParameter');
    switch currentNode.schedStrategy
        case SchedStrategy.PS
            subParameterNode.setAttribute('classPath', 'jmt.engine.NetStrategies.PSStrategies.EPSStrategy');
            subParameterNode.setAttribute('name', 'EPSStrategy');
        case SchedStrategy.DPS
            subParameterNode.setAttribute('classPath', 'jmt.engine.NetStrategies.PSStrategies.DPSStrategy');
            subParameterNode.setAttribute('name', 'DPSStrategy');
        case SchedStrategy.GPS
            subParameterNode.setAttribute('classPath', 'jmt.engine.NetStrategies.PSStrategies.GPSStrategy');
            subParameterNode.setAttribute('name', 'GPSStrategy');
    end
    
    visitsNode.appendChild(subParameterNode);
    section.appendChild(visitsNode);
end
end

