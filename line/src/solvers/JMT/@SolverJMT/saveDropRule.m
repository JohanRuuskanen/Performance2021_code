function [simDoc, section] = saveDropRule(self, simDoc, section, currentNode)
    % [SIMDOC, SECTION] = SAVEDROPRULE(SIMDOC, SECTION)
    
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    schedStrategyNode = simDoc.createElement('parameter');
    schedStrategyNode.setAttribute('array', 'true');
    schedStrategyNode.setAttribute('classPath', 'java.lang.String');
    schedStrategyNode.setAttribute('name', 'dropRules');
    numOfClasses = length(self.model.classes);
    for i=1:(numOfClasses)
        currentClass = self.model.classes{i,1};
        
        refClassNode = simDoc.createElement('refClass');
        refClassNode.appendChild(simDoc.createTextNode(currentClass.name));
        schedStrategyNode.appendChild(refClassNode);
        
        subParameterNode = simDoc.createElement('subParameter');
        subParameterNode.setAttribute('classPath', 'java.lang.String');
        subParameterNode.setAttribute('name', 'dropRule');
        
        valueNode2 = simDoc.createElement('value');
        valueNode2.appendChild(simDoc.createTextNode(DropStrategy.toText(currentNode.droprule(i))));
        
        subParameterNode.appendChild(valueNode2);
        schedStrategyNode.appendChild(subParameterNode);
        section.appendChild(schedStrategyNode);
    end
    end
    