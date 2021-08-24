function [simDoc, section] = saveClassSwitchStrategy(self, simDoc, section, currentNode)
% [SIMDOC, SECTION] = SAVECLASSSWITCHSTRATEGY(SIMDOC, SECTION, CURRENTNODE)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

paramNode = simDoc.createElement('parameter');
paramNode.setAttribute('array', 'true');
paramNode.setAttribute('classPath', 'java.lang.Object');
paramNode.setAttribute('name', 'matrix');

numOfClasses = length(self.model.classes);
for i=1:numOfClasses
    currentClass = self.model.classes{i,1};
    
    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(currentClass.name));
    paramNode.appendChild(refClassNode);
    
    
    subParNodeRow = simDoc.createElement('subParameter');
    subParNodeRow.setAttribute('array', 'true');
    subParNodeRow.setAttribute('classPath', 'java.lang.Float');
    subParNodeRow.setAttribute('name', 'row');
    for j=1:numOfClasses
        nextClass = self.model.classes{j,1};
        
        refClassNode = simDoc.createElement('refClass');
        refClassNode.appendChild(simDoc.createTextNode(nextClass.name));
        subParNodeRow.appendChild(refClassNode);
        
        subParNodeCell = simDoc.createElement('subParameter');
        subParNodeCell.setAttribute('classPath', 'java.lang.Float');
        subParNodeCell.setAttribute('name', 'cell');
        valNode = simDoc.createElement('value');
        valNode.appendChild(simDoc.createTextNode(sprintf('%12.12f',currentNode.server.csFun(i,j))));
        subParNodeCell.appendChild(valNode);
        subParNodeRow.appendChild(subParNodeCell);
        
    end
    paramNode.appendChild(subParNodeRow);
    
end
section.appendChild(paramNode);

end
