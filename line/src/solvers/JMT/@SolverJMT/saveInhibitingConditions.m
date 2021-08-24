function [simDoc, section] = saveInhibitingConditions(self, simDoc, section, currentNode)
    % [SIMDOC, SECTION] = SAVEFIRINGOUTCOMES(SIMDOC, SECTION, CURRENTNODE)
    
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    inhibitingConditionsNode = simDoc.createElement('parameter');
    inhibitingConditionsNode.setAttribute('array', 'true');
    inhibitingConditionsNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix');
    inhibitingConditionsNode.setAttribute('name', 'inhibitingConditions');

    numOflinks = length(self.model.links);
    connections = {};
    inputs = [];
    for j=1:(numOflinks)
        currentConnection = self.model.links{j,1};
        if strcmp(currentConnection{2}.name, currentNode.name)
            inputName = currentConnection{1}.name;
            connections{end+1} = inputName;
            inputs = [inputs, self.model.getNodeIndex(inputName)];
        end
    end
    numOfInputs = length(connections);

    numOfModes = length(currentNode.modeNames);
    numOfClasses = length(self.model.classes);
    for i=1:(numOfModes)

        subInhibitingConditionNode = simDoc.createElement('subParameter');
        subInhibitingConditionNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix');
        subInhibitingConditionNode.setAttribute('name', 'inhibitingCondition');

        subInhibitingVectorsNode = simDoc.createElement('subParameter');
        subInhibitingVectorsNode.setAttribute('array', 'true');
        subInhibitingVectorsNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector');
        subInhibitingVectorsNode.setAttribute('name', 'inhibitingVectors');

        for k=1:(numOfInputs)
            subInhibitingVectorNode = simDoc.createElement('subParameter');
            subInhibitingVectorNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector');
            subInhibitingVectorNode.setAttribute('name', 'inhibitingVector');

            subStationNameNode = simDoc.createElement('subParameter');
            subStationNameNode.setAttribute('classPath', 'java.lang.String');
            subStationNameNode.setAttribute('name', 'stationName');

            placeNameValueNode = simDoc.createElement('value');
            placeNameValueNode.appendChild(simDoc.createTextNode(connections(k)));
            subStationNameNode.appendChild(placeNameValueNode);

            subInhibitingVectorNode.appendChild(subStationNameNode);

            subInhibitingEntriesNode = simDoc.createElement('subParameter');
            subInhibitingEntriesNode.setAttribute('array', 'true');
            subInhibitingEntriesNode.setAttribute('classPath', 'java.lang.Integer');
            subInhibitingEntriesNode.setAttribute('name', 'inhibitingEntries');
            
            for j=1:(numOfClasses)
                currentClass = self.model.classes{j,1};
            
                refClassNode = simDoc.createElement('refClass');
                refClassNode.appendChild(simDoc.createTextNode(currentClass.name));
                subInhibitingEntriesNode.appendChild(refClassNode);
    
                subParameterNode = simDoc.createElement('subParameter');
                subParameterNode.setAttribute('classPath', 'java.lang.Integer');
                subParameterNode.setAttribute('name', 'inhibitingEntry');
                
                valueNode2 = simDoc.createElement('value');

                if isinf(currentNode.inhibitingConditions(k,j,i))
                    valueNode2.appendChild(simDoc.createTextNode(int2str(-1)));
                else
                    valueNode2.appendChild(simDoc.createTextNode(int2str(currentNode.inhibitingConditions(inputs(k),j,i))));
                end

                subParameterNode.appendChild(valueNode2);
                subInhibitingEntriesNode.appendChild(subParameterNode);
                subInhibitingVectorNode.appendChild(subInhibitingEntriesNode);
            end
            subInhibitingVectorsNode.appendChild(subInhibitingVectorNode);
        end
        
        subInhibitingConditionNode.appendChild(subInhibitingVectorsNode);
        inhibitingConditionsNode.appendChild(subInhibitingConditionNode);
    end
    section.appendChild(inhibitingConditionsNode);
    end