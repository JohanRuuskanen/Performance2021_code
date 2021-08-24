function [simDoc, section] = saveEnablingConditions(self, simDoc, section, currentNode)
    % [SIMDOC, SECTION] = SAVEFIRINGOUTCOMES(SIMDOC, SECTION, CURRENTNODE)
    
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    enablingConditionsNode = simDoc.createElement('parameter');
    enablingConditionsNode.setAttribute('array', 'true');
    enablingConditionsNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix');
    enablingConditionsNode.setAttribute('name', 'enablingConditions');

    % numOflinks = length(self.model.links);
    % connections = {};
    % inputs = [];
    % for j=1:(numOflinks)
    %     currentConnection = self.model.links{j,1};
    %     if strcmp(currentConnection{2}.name, currentNode.name)
    %         inputName = currentConnection{1}.name;
    %         connections{end+1} = inputName;
    %         inputs = [inputs, self.model.getNodeIndex(inputName)];
    %     end
    % end
    % numOfInputs = length(connections);
    numOfNodes = length(self.model.nodes);
    numOfModes = length(currentNode.modeNames);
    numOfClasses = length(self.model.classes);
    for i=1:(numOfModes)

        subEnablingConditionNode = simDoc.createElement('subParameter');
        subEnablingConditionNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix');
        subEnablingConditionNode.setAttribute('name', 'enablingCondition');

        subEnablingVectorsNode = simDoc.createElement('subParameter');
        subEnablingVectorsNode.setAttribute('array', 'true');
        subEnablingVectorsNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector');
        subEnablingVectorsNode.setAttribute('name', 'enablingVectors');

        for k=1:(numOfNodes)
            subEnablingVectorNode = simDoc.createElement('subParameter');
            subEnablingVectorNode.setAttribute('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector');
            subEnablingVectorNode.setAttribute('name', 'enablingVector');

            subStationNameNode = simDoc.createElement('subParameter');
            subStationNameNode.setAttribute('classPath', 'java.lang.String');
            subStationNameNode.setAttribute('name', 'stationName');

            placeNameValueNode = simDoc.createElement('value');
            placeNameValueNode.appendChild(simDoc.createTextNode(self.model.getNodeByIndex(k).name));
            subStationNameNode.appendChild(placeNameValueNode);

            subEnablingVectorNode.appendChild(subStationNameNode);

            subEnablingEntriesNode = simDoc.createElement('subParameter');
            subEnablingEntriesNode.setAttribute('array', 'true');
            subEnablingEntriesNode.setAttribute('classPath', 'java.lang.Integer');
            subEnablingEntriesNode.setAttribute('name', 'enablingEntries');
    
            exists = false;

            for j=1:(numOfClasses)
                currentClass = self.model.classes{j,1};
                refClassNode = simDoc.createElement('refClass');
                refClassNode.appendChild(simDoc.createTextNode(currentClass.name));
                subEnablingEntriesNode.appendChild(refClassNode);
    
                subParameterNode = simDoc.createElement('subParameter');
                subParameterNode.setAttribute('classPath', 'java.lang.Integer');
                subParameterNode.setAttribute('name', 'enablingEntry');
                
                valueNode2 = simDoc.createElement('value');

                if isinf(currentNode.enablingConditions(k,j,i))
                    valueNode2.appendChild(simDoc.createTextNode(int2str(-1)));
                    exists = true;
                elseif currentNode.enablingConditions(k,j,i) > 0
                    valueNode2.appendChild(simDoc.createTextNode(int2str(currentNode.enablingConditions(k,j,i))));
                    exists = true;
                elseif ~isinf(currentNode.inhibitingConditions(k,j,i)) && currentNode.inhibitingConditions(k,j,i) > 0
                    valueNode2.appendChild(simDoc.createTextNode(int2str(0)));
                    exists = true;
                end

                subParameterNode.appendChild(valueNode2);
                subEnablingEntriesNode.appendChild(subParameterNode);
                subEnablingVectorNode.appendChild(subEnablingEntriesNode);
            end
            if exists
                subEnablingVectorsNode.appendChild(subEnablingVectorNode);
            end
            
        end
        
        subEnablingConditionNode.appendChild(subEnablingVectorsNode);
        enablingConditionsNode.appendChild(subEnablingConditionNode);
    end
    section.appendChild(enablingConditionsNode);
    end