function [simDoc, section] = savePreemptiveWeights(self, simDoc, section, currentNode)
% [SIMDOC, SECTION] = SAVEPREEMPTIVEWEIGHTS(SIMDOC, SECTION, CURRENTNODE)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
visitsNode = simDoc.createElement('parameter');
visitsNode.setAttribute('array', 'true');
visitsNode.setAttribute('classPath', 'java.lang.Double');
visitsNode.setAttribute('name', 'serviceWeights');

numOfClasses = length(self.model.classes);
for i=1:(numOfClasses)
    currentClass = self.model.classes{i,1};
    
    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(currentClass.name));
    visitsNode.appendChild(refClassNode);
    
    subParameterNode = simDoc.createElement('subParameter');
    subParameterNode.setAttribute('classPath', 'java.lang.Double');
    subParameterNode.setAttribute('name', 'serviceWeight');
    
    valueNode2 = simDoc.createElement('value');
    if isempty(currentNode.schedStrategyPar) % PS case
        valueNode2.appendChild(simDoc.createTextNode(int2str(1)));
    else
        valueNode2.appendChild(simDoc.createTextNode(num2str(currentNode.schedStrategyPar(i))));
    end
    
    subParameterNode.appendChild(valueNode2);
    visitsNode.appendChild(subParameterNode);
    section.appendChild(visitsNode);
end
end
