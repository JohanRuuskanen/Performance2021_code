function [simDoc, section] = saveDropStrategy(self, simDoc, section)
% [SIMDOC, SECTION] = SAVEDROPSTRATEGY(SIMDOC, SECTION)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

schedStrategyNode = simDoc.createElement('parameter');
schedStrategyNode.setAttribute('array', 'true');
schedStrategyNode.setAttribute('classPath', 'java.lang.String');
schedStrategyNode.setAttribute('name', 'dropStrategies');
numOfClasses = length(self.model.classes);
for i=1:(numOfClasses)
    currentClass = self.model.classes{i,1};
    
    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(currentClass.name));
    schedStrategyNode.appendChild(refClassNode);
    
    subParameterNode = simDoc.createElement('subParameter');
    subParameterNode.setAttribute('classPath', 'java.lang.String');
    subParameterNode.setAttribute('name', 'dropStrategy');
    
    valueNode2 = simDoc.createElement('value');
    valueNode2.appendChild(simDoc.createTextNode('drop'));
    
    subParameterNode.appendChild(valueNode2);
    schedStrategyNode.appendChild(subParameterNode);
    section.appendChild(schedStrategyNode);
end
end
