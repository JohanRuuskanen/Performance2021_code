function [simDoc, section] = savePlaceCapacities(self, simDoc, section, currentNode)
    % [SIMDOC, SECTION] = SAVEPLACECAPACITY(SIMDOC, SECTION)
    
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    placeCapacityNode = simDoc.createElement('parameter');
    placeCapacityNode.setAttribute('array', 'true');
    placeCapacityNode.setAttribute('classPath', 'java.lang.Integer');
    placeCapacityNode.setAttribute('name', 'capacities');
    numOfClasses = length(self.model.classes);
    for i=1:(numOfClasses)
        currentClass = self.model.classes{i,1};
        
        refClassNode = simDoc.createElement('refClass');
        refClassNode.appendChild(simDoc.createTextNode(currentClass.name));
        placeCapacityNode.appendChild(refClassNode);
        
        subParameterNode = simDoc.createElement('subParameter');
        subParameterNode.setAttribute('classPath', 'java.lang.Integer');
        subParameterNode.setAttribute('name', 'capacity');
        
        valueNode2 = simDoc.createElement('value');
        if isinf(currentNode.classCap(i))
            valueNode2.appendChild(simDoc.createTextNode(int2str(-1)));
        else
            valueNode2.appendChild(simDoc.createTextNode(int2str(currentNode.classCap(i))));
        end
        
        subParameterNode.appendChild(valueNode2);
        placeCapacityNode.appendChild(subParameterNode);
    end
    section.appendChild(placeCapacityNode);
    end
    