function [simDoc, section] = saveNumbersOfServers(self, simDoc, section, currentNode)
    % [SIMDOC, SECTION] = SAVENUMBERSOFSERVERS(SIMDOC, SECTION, CURRENTNODE)
    
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.

    numbersOfServersNode = simDoc.createElement('parameter');
    numbersOfServersNode.setAttribute('classPath', 'java.lang.Integer');
    numbersOfServersNode.setAttribute('name', 'numbersOfServers');
    numbersOfServersNode.setAttribute('array', 'true');

    numOfModes = length(currentNode.modeNames);
    for i=1:(numOfModes)

        subNumberOfServersNode = simDoc.createElement('subParameter');
        subNumberOfServersNode.setAttribute('classPath', 'java.lang.Integer');
        subNumberOfServersNode.setAttribute('name', 'numberOfServers');
    
        valueNode = simDoc.createElement('value');

        if isinf(currentNode.numbersOfServers(i))
            valueNode.appendChild(simDoc.createTextNode(int2str(-1)));
        else
            valueNode.appendChild(simDoc.createTextNode(int2str(currentNode.numbersOfServers(i))));
        end

        subNumberOfServersNode.appendChild(valueNode);
        numbersOfServersNode.appendChild(subNumberOfServersNode);
    end

    section.appendChild(numbersOfServersNode);
    end
    