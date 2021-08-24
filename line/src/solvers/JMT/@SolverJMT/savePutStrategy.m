function [simDoc, section] = savePutStrategy(self, simDoc, section, currentNode)
% [SIMDOC, SECTION] = SAVEPUTSTRATEGY(SIMDOC, SECTION, CURRENTNODE)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
queuePutStrategyNode = simDoc.createElement('parameter');
queuePutStrategyNode.setAttribute('array', 'true');
queuePutStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategy');
queuePutStrategyNode.setAttribute('name', 'QueuePutStrategy');

numOfClasses = length(self.model.classes);
for i=1:(numOfClasses)
    currentClass = self.model.classes{i,1};
    
    refClassNode2 = simDoc.createElement('refClass');
    refClassNode2.appendChild(simDoc.createTextNode(currentClass.name));
    
    queuePutStrategyNode.appendChild(refClassNode2);
    %                switch currentNode.input.inputJobClasses{i}{2}
    switch currentNode.schedStrategy
        case SchedStrategy.SIRO
            subParameterNode2 = simDoc.createElement('subParameter');
            subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.RandStrategy');
            subParameterNode2.setAttribute('name', 'RandStrategy');
        case SchedStrategy.LJF
            subParameterNode2 = simDoc.createElement('subParameter');
            subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.LJFStrategy');
            subParameterNode2.setAttribute('name', 'LJFStrategy');
        case SchedStrategy.SJF
            subParameterNode2 = simDoc.createElement('subParameter');
            subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.SJFStrategy');
            subParameterNode2.setAttribute('name', 'SJFStrategy');
        case SchedStrategy.LEPT
            subParameterNode2 = simDoc.createElement('subParameter');
            subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.LEPTStrategy');
            subParameterNode2.setAttribute('name', 'LEPTStrategy');
        case SchedStrategy.SEPT
            subParameterNode2 = simDoc.createElement('subParameter');
            subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.SEPTStrategy');
            subParameterNode2.setAttribute('name', 'SEPTStrategy');
        case SchedStrategy.LCFS
            subParameterNode2 = simDoc.createElement('subParameter');
            subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.HeadStrategy');
            subParameterNode2.setAttribute('name', 'HeadStrategy');
        case SchedStrategy.HOL
            subParameterNode2 = simDoc.createElement('subParameter');
            subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategyPriority');
            subParameterNode2.setAttribute('name', 'TailStrategy');
        otherwise % treat as FCFS - this is required for PS
            subParameterNode2 = simDoc.createElement('subParameter');
            subParameterNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy');
            subParameterNode2.setAttribute('name', 'TailStrategy');
    end
    queuePutStrategyNode.appendChild(subParameterNode2);
    section.appendChild(queuePutStrategyNode);
end
end
