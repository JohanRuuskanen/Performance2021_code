function [simDoc, section] = saveServerVisits(self, simDoc, section)
% [SIMDOC, SECTION] = SAVESERVERVISITS(SIMDOC, SECTION)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
visitsNode = simDoc.createElement('parameter');
visitsNode.setAttribute('array', 'true');
visitsNode.setAttribute('classPath', 'java.lang.Integer');
visitsNode.setAttribute('name', 'numberOfVisits');

numOfClasses = length(self.model.classes);
for i=1:(numOfClasses)
    currentClass = self.model.classes{i,1};
    
    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(currentClass.name));
    visitsNode.appendChild(refClassNode);
    
    subParameterNode = simDoc.createElement('subParameter');
    subParameterNode.setAttribute('classPath', 'java.lang.Integer');
    subParameterNode.setAttribute('name', 'numberOfVisits');
    
    valueNode2 = simDoc.createElement('value');
    valueNode2.appendChild(simDoc.createTextNode(int2str(1)));
    
    subParameterNode.appendChild(valueNode2);
    visitsNode.appendChild(subParameterNode);
    section.appendChild(visitsNode);
end
end
