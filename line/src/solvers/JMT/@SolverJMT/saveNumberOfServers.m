function [simDoc, section] = saveNumberOfServers(self, simDoc, section, currentNode)
% [SIMDOC, SECTION] = SAVENUMBEROFSERVERS(SIMDOC, SECTION, CURRENTNODE)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
sizeNode = simDoc.createElement('parameter');
sizeNode.setAttribute('classPath', 'java.lang.Integer');
sizeNode.setAttribute('name', 'maxJobs');

valueNode = simDoc.createElement('value');
valueNode.appendChild(simDoc.createTextNode(int2str(currentNode.numberOfServers)));

sizeNode.appendChild(valueNode);
section.appendChild(sizeNode);
end
