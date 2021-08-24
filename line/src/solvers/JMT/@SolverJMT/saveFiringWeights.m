function [simDoc, section] = saveFiringWeights(self, simDoc, section, currentNode)
    % [SIMDOC, SECTION] = SAVEFIRINGWEIGHTS(SIMDOC, SECTION, CURRENTNODE)
    
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.

    firingWeightsNode = simDoc.createElement('parameter');
    firingWeightsNode.setAttribute('classPath', 'java.lang.Double');
    firingWeightsNode.setAttribute('name', 'firingWeights');
    firingWeightsNode.setAttribute('array', 'true');

    numOfModes = length(currentNode.modeNames);
    for i=1:(numOfModes)

        subFiringWeightNode = simDoc.createElement('subParameter');
        subFiringWeightNode.setAttribute('classPath', 'java.lang.Double');
        subFiringWeightNode.setAttribute('name', 'firingWeight');
    
        valueNode = simDoc.createElement('value');
        valueNode.appendChild(simDoc.createTextNode(int2str(currentNode.firingWeights(i))));

        subFiringWeightNode.appendChild(valueNode);
        firingWeightsNode.appendChild(subFiringWeightNode);
    end

    section.appendChild(firingWeightsNode);
    end
    