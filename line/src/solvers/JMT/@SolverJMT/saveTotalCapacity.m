function [simDoc, section] = saveTotalCapacity(self, simDoc, section, currentNode)
    % [SIMDOC, SECTION] = SAVETOTALCAPACITY(SIMDOC, SECTION, CURRENTNODE)
    
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    sizeNode = simDoc.createElement('parameter');
    sizeNode.setAttribute('classPath', 'java.lang.Integer');
    sizeNode.setAttribute('name', 'totalCapacity');
    valueNode = simDoc.createElement('value');
    if isinf(currentNode.cap)
        valueNode.appendChild(simDoc.createTextNode(int2str(-1)));
    else
        valueNode.appendChild(simDoc.createTextNode(int2str(currentNode.cap)));
    end
    
    sizeNode.appendChild(valueNode);
    section.appendChild(sizeNode);
    end
    