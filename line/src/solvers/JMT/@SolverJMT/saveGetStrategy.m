function [simDoc, section] = saveGetStrategy(self, simDoc, section, currentNode)
% [SIMDOC, SECTION] = SAVEGETSTRATEGY(SIMDOC, SECTION, CURRENTNODE)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

% the get strategy is always fcfs
queueGetStrategyNode = simDoc.createElement('parameter');
queueGetStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy');
queueGetStrategyNode.setAttribute('name', 'FCFSstrategy');
section.appendChild(queueGetStrategyNode);
end
