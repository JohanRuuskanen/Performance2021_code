function [simDoc, section] = saveRoutingStrategy(self, simDoc, section, currentNode)
% [SIMDOC, SECTION] = SAVEROUTINGSTRATEGY(SIMDOC, SECTION, CURRENTNODE)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
strategyNode = simDoc.createElement('parameter');
strategyNode.setAttribute('array', 'true');
strategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.RoutingStrategy');
strategyNode.setAttribute('name', 'RoutingStrategy');

numOfClasses = length(self.model.classes);
for i=1:(numOfClasses)
    currentClass = self.model.classes{i,1};
    
    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(currentClass.name));
    strategyNode.appendChild(refClassNode);
    
    switch currentNode.output.outputStrategy{i}{2}
        case RoutingStrategy.RAND
            concStratNode = simDoc.createElement('subParameter');
            concStratNode.setAttribute('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.RandomStrategy');
            concStratNode.setAttribute('name', 'Random');
        case RoutingStrategy.RRB
            concStratNode = simDoc.createElement('subParameter');
            concStratNode.setAttribute('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.RoundRobinStrategy');
            concStratNode.setAttribute('name', 'Round Robin');
        case RoutingStrategy.JSQ
            concStratNode = simDoc.createElement('subParameter');
            concStratNode.setAttribute('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.ShortestQueueLengthRoutingStrategy');
            concStratNode.setAttribute('name', 'Join the Shortest Queue (JSQ)');
        case RoutingStrategy.PROB
            concStratNode = simDoc.createElement('subParameter');
            concStratNode.setAttribute('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.EmpiricalStrategy');
            concStratNode.setAttribute('name', RoutingStrategy.PROB);
            concStratNode2 = simDoc.createElement('subParameter');
            concStratNode2.setAttribute('array', 'true');
            concStratNode2.setAttribute('classPath', 'jmt.engine.random.EmpiricalEntry');
            concStratNode2.setAttribute('name', 'EmpiricalEntryArray');
            for k=1:length(currentNode.output.outputStrategy{i}{end})
                concStratNode3 = simDoc.createElement('subParameter');
                concStratNode3.setAttribute('classPath', 'jmt.engine.random.EmpiricalEntry');
                concStratNode3.setAttribute('name', 'EmpiricalEntry');
                concStratNode4Station = simDoc.createElement('subParameter');
                concStratNode4Station.setAttribute('classPath', 'java.lang.String');
                concStratNode4Station.setAttribute('name', 'stationName');
                concStratNode4StationValueNode = simDoc.createElement('value');
                concStratNode4StationValueNode.appendChild(simDoc.createTextNode(sprintf('%s',currentNode.output.outputStrategy{i}{end}{k}{1}.name)));
                concStratNode4Station.appendChild(concStratNode4StationValueNode);
                concStratNode3.appendChild(concStratNode4Station);
                concStratNode4Probability = simDoc.createElement('subParameter');
                concStratNode4Probability.setAttribute('classPath', 'java.lang.Double');
                concStratNode4Probability.setAttribute('name', 'probability');
                concStratNode4ProbabilityValueNode = simDoc.createElement('value');
                concStratNode4ProbabilityValueNode.appendChild(simDoc.createTextNode(sprintf('%12.12f',currentNode.output.outputStrategy{i}{end}{k}{2})));
                concStratNode4Probability.appendChild(concStratNode4ProbabilityValueNode);
                
                concStratNode3.appendChild(concStratNode4Station);
                concStratNode3.appendChild(concStratNode4Probability);
                concStratNode2.appendChild(concStratNode3);
            end
            concStratNode.appendChild(concStratNode2);
        otherwise
            concStratNode = simDoc.createElement('subParameter');
            concStratNode.setAttribute('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.DisabledRoutingStrategy');
            concStratNode.setAttribute('name', 'Random');
    end
    strategyNode.appendChild(concStratNode);
    section.appendChild(strategyNode);
end
end
