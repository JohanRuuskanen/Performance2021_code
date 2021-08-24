function [simDoc, section] = saveForkStrategy(self, simDoc, section, currentNode)
% [SIMDOC, SECTION] = SAVEFORKSTRATEGY(SIMDOC, SECTION, CURRENTNODE)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

jplNode = simDoc.createElement('parameter');
jplNode.setAttribute('classPath', 'java.lang.Integer');
jplNode.setAttribute('name', 'jobsPerLink');
valueNode = simDoc.createElement('value');
valueNode.appendChild(simDoc.createTextNode(int2str(currentNode.output.tasksPerLink)));
jplNode.appendChild(valueNode);
section.appendChild(jplNode);

blockNode = simDoc.createElement('parameter');
blockNode.setAttribute('classPath', 'java.lang.Integer');
blockNode.setAttribute('name', 'block');
valueNode = simDoc.createElement('value');
valueNode.appendChild(simDoc.createTextNode(int2str(-1)));
blockNode.appendChild(valueNode);
section.appendChild(blockNode);

issimplNode = simDoc.createElement('parameter');
issimplNode.setAttribute('classPath', 'java.lang.Boolean');
issimplNode.setAttribute('name', 'isSimplifiedFork');
valueNode = simDoc.createElement('value');
valueNode.appendChild(simDoc.createTextNode('true'));
issimplNode.appendChild(valueNode);
section.appendChild(issimplNode);

strategyNode = simDoc.createElement('parameter');
strategyNode.setAttribute('array', 'true');
strategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ForkStrategy');
strategyNode.setAttribute('name', 'ForkStrategy');

numOfClasses = length(self.model.classes);
for class=1:(numOfClasses)
    currentClass = self.model.classes{class,1};
    
    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(currentClass.name));
    strategyNode.appendChild(refClassNode);
    
    classStratNode = simDoc.createElement('subParameter');
    classStratNode.setAttribute('classPath', 'jmt.engine.NetStrategies.ForkStrategies.ProbabilitiesFork');
    classStratNode.setAttribute('name', 'Branch Probabilities');
    classStratNode2 = simDoc.createElement('subParameter');
    classStratNode2.setAttribute('array', 'true');
    classStratNode2.setAttribute('classPath', 'jmt.engine.NetStrategies.ForkStrategies.OutPath');
    classStratNode2.setAttribute('name', 'EmpiricalEntryArray');
    switch currentNode.output.outputStrategy{class}{2} % check routing strategy
        case RoutingStrategy.PROB
            for k=1:length(currentNode.output.outputStrategy{class}{3})
                classStratNode3 = simDoc.createElement('subParameter');
                classStratNode3.setAttribute('classPath', 'jmt.engine.NetStrategies.ForkStrategies.OutPath');
                classStratNode3.setAttribute('name', 'OutPathEntry');
                classStratNode4 = simDoc.createElement('subParameter');
                classStratNode4.setAttribute('classPath', 'jmt.engine.random.EmpiricalEntry');
                classStratNode4.setAttribute('name', 'outUnitProbability');
                classStratNode4Station = simDoc.createElement('subParameter');
                classStratNode4Station.setAttribute('classPath', 'java.lang.String');
                classStratNode4Station.setAttribute('name', 'stationName');
                classStratNode4StationValueNode = simDoc.createElement('value');
                %            classStratNode4StationValueNode.appendChild(simDoc.createTextNode(sprintf('%s',RoutingStrategy.toText(currentNode.output.outputStrategy{i}{end}{k}{2}))));
                %            classStratNode4StationValueNode.appendChild(simDoc.createTextNode(sprintf('%s',RoutingStrategy.toText(currentNode.output.outputStrategy{i}{2}))));
                classStratNode4StationValueNode.appendChild(simDoc.createTextNode(sprintf('%s',currentNode.output.outputStrategy{class}{end}{k}{1}.name)));
            end
            classStratNode4Station.appendChild(classStratNode4StationValueNode);
            classStratNode3.appendChild(classStratNode4Station);
            classStratNode4Probability = simDoc.createElement('subParameter');
            classStratNode4Probability.setAttribute('classPath', 'java.lang.Double');
            classStratNode4Probability.setAttribute('name', 'probability');
            classStratNode4ProbabilityValueNode = simDoc.createElement('value');
            classStratNode4ProbabilityValueNode.appendChild(simDoc.createTextNode('1.0'));
            classStratNode4Probability.appendChild(classStratNode4ProbabilityValueNode);
            
            classStratNode4b = simDoc.createElement('subParameter');
            classStratNode4b.setAttribute('classPath', 'jmt.engine.random.EmpiricalEntry');
            classStratNode4b.setAttribute('array', 'true');
            classStratNode4b.setAttribute('name', 'JobsPerLinkDis');
            classStratNode5b = simDoc.createElement('subParameter');
            classStratNode5b.setAttribute('classPath', 'jmt.engine.random.EmpiricalEntry');
            classStratNode5b.setAttribute('name', 'EmpiricalEntry');
            classStratNode5bStation = simDoc.createElement('subParameter');
            classStratNode5bStation.setAttribute('classPath', 'java.lang.String');
            classStratNode5bStation.setAttribute('name', 'numbers');
            classStratNode5bStationValueNode = simDoc.createElement('value');
            classStratNode5bStationValueNode.appendChild((simDoc.createTextNode(int2str(currentNode.output.tasksPerLink))));
            classStratNode5bStation.appendChild(classStratNode5bStationValueNode);
            classStratNode4b.appendChild(classStratNode5bStation);
            classStratNode5bProbability = simDoc.createElement('subParameter');
            classStratNode5bProbability.setAttribute('classPath', 'java.lang.Double');
            classStratNode5bProbability.setAttribute('name', 'probability');
            classStratNode5bProbabilityValueNode = simDoc.createElement('value');
            classStratNode5bProbabilityValueNode.appendChild(simDoc.createTextNode('1.0'));
            classStratNode5bProbability.appendChild(classStratNode5bProbabilityValueNode);
            
            classStratNode4.appendChild(classStratNode4Station);
            classStratNode4.appendChild(classStratNode4Probability);
            classStratNode3.appendChild(classStratNode4);
            classStratNode5b.appendChild(classStratNode5bStation);
            classStratNode5b.appendChild(classStratNode5bProbability);
            classStratNode4b.appendChild(classStratNode5b);
            classStratNode3.appendChild(classStratNode4b);
            classStratNode2.appendChild(classStratNode3);
    end
    classStratNode.appendChild(classStratNode2);
    strategyNode.appendChild(classStratNode);
end
section.appendChild(strategyNode);
end
