function [simDoc, section] = saveFiringPriorities(self, simDoc, section, currentNode)
    % [SIMDOC, SECTION] = SAVEFIRINGPRIORITIES(SIMDOC, SECTION, CURRENTNODE)
    
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.

    firingPrioritiesNode = simDoc.createElement('parameter');
    firingPrioritiesNode.setAttribute('classPath', 'java.lang.Integer');
    firingPrioritiesNode.setAttribute('name', 'firingPriorities');
    firingPrioritiesNode.setAttribute('array', 'true');

    numOfModes = length(currentNode.modeNames);
    for i=1:(numOfModes)

        subFiringPriorityNode = simDoc.createElement('subParameter');
        subFiringPriorityNode.setAttribute('classPath', 'java.lang.Integer');
        subFiringPriorityNode.setAttribute('name', 'firingPriority');
    
        valueNode = simDoc.createElement('value');
        if isinf(currentNode.firingPriorities(i))
            valueNode.appendChild(simDoc.createTextNode(int2str(-1)));
        else
            valueNode.appendChild(simDoc.createTextNode(int2str(currentNode.firingPriorities(i))));
        end
        
        subFiringPriorityNode.appendChild(valueNode);
        firingPrioritiesNode.appendChild(subFiringPriorityNode);
    end

    section.appendChild(firingPrioritiesNode);
    end
    