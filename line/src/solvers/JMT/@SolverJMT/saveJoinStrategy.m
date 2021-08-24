function [simDoc, section] = saveJoinStrategy(self, simDoc, section, currentNode)
% [SIMDOC, SECTION] = SAVEJOINSTRATEGY(SIMDOC, SECTION, CURRENTNODE)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
strategyNode = simDoc.createElement('parameter');
strategyNode.setAttribute('array', 'true');
strategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.JoinStrategy');
strategyNode.setAttribute('name', 'JoinStrategy');

numOfClasses = length(self.model.classes);
for i=1:(numOfClasses)
    currentClass = self.model.classes{i,1};
    switch currentNode.input.joinStrategy{currentClass.index}
        case JoinStrategy.STD
            refClassNode2 = simDoc.createElement('refClass');
            refClassNode2.appendChild(simDoc.createTextNode(currentClass.name));
            strategyNode.appendChild(refClassNode2);
            
            joinStrategyNode = simDoc.createElement('subParameter');
            joinStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.JoinStrategies.NormalJoin');
            joinStrategyNode.setAttribute('name', 'Standard Join');
            reqNode = simDoc.createElement('subParameter');
            reqNode.setAttribute('classPath', 'java.lang.Integer');
            reqNode.setAttribute('name', 'numRequired');
            valueNode = simDoc.createElement('value');
            valueNode.appendChild(simDoc.createTextNode(int2str(currentNode.input.joinRequired{currentClass.index})));
            reqNode.appendChild(valueNode);
            joinStrategyNode.appendChild(reqNode);
            strategyNode.appendChild(joinStrategyNode);
            section.appendChild(strategyNode);
        case JoinStrategy.Quorum
            refClassNode2 = simDoc.createElement('refClass');
            refClassNode2.appendChild(simDoc.createTextNode(currentClass.name));
            strategyNode.appendChild(refClassNode2);
            
            joinStrategyNode = simDoc.createElement('subParameter');
            joinStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.JoinStrategies.PartialJoin');
            joinStrategyNode.setAttribute('name', 'Quorum');
            reqNode = simDoc.createElement('subParameter');
            reqNode.setAttribute('classPath', 'java.lang.Integer');
            reqNode.setAttribute('name', 'numRequired');
            valueNode = simDoc.createElement('value');
            valueNode.appendChild(simDoc.createTextNode(int2str(currentNode.input.joinRequired{currentClass.index})));
            reqNode.appendChild(valueNode);
            joinStrategyNode.appendChild(reqNode);
            strategyNode.appendChild(joinStrategyNode);
            section.appendChild(strategyNode);
        case JoinStrategy.Guard
            refClassNode2 = simDoc.createElement('refClass');
            refClassNode2.appendChild(simDoc.createTextNode(currentClass.name));
            strategyNode.appendChild(refClassNode2);
            
            joinStrategyNode = simDoc.createElement('subParameter');
            joinStrategyNode.setAttribute('classPath', 'jmt.engine.NetStrategies.JoinStrategies.PartialJoin');
            joinStrategyNode.setAttribute('name', 'Quorum');
            reqNode = simDoc.createElement('subParameter');
            reqNode.setAttribute('classPath', 'java.lang.Integer');
            reqNode.setAttribute('name', 'numRequired');
            valueNode = simDoc.createElement('value');
            valueNode.appendChild(simDoc.createTextNode(int2str(currentNode.input.joinRequired{currentClass.index})));
            reqNode.appendChild(valueNode);
            joinStrategyNode.appendChild(reqNode);
            strategyNode.appendChild(joinStrategyNode);
            section.appendChild(strategyNode);
    end
end
end
