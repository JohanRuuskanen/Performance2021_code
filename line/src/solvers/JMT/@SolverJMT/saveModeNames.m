function [simDoc, section] = saveModeNames(self, simDoc, section, currentNode)
    % [SIMDOC, SECTION] = SAVEMODENAMES(SIMDOC, SECTION, CURRENTNODE)
    
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.

    modeNamesNode = simDoc.createElement('parameter');
    modeNamesNode.setAttribute('classPath', 'java.lang.String');
    modeNamesNode.setAttribute('name', 'modeNames');
    modeNamesNode.setAttribute('array', 'true');

    numOfModes = length(currentNode.modeNames);
    for i=1:(numOfModes)

        subModeNameNode = simDoc.createElement('subParameter');
        subModeNameNode.setAttribute('classPath', 'java.lang.String');
        subModeNameNode.setAttribute('name', 'modeName');
    
        valueNode = simDoc.createElement('value');
        valueNode.appendChild(simDoc.createTextNode(currentNode.modeNames{i}));

        subModeNameNode.appendChild(valueNode);
        modeNamesNode.appendChild(subModeNameNode);
    end

    section.appendChild(modeNamesNode);
    end
    