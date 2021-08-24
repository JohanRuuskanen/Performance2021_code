function fname = writeJSIM(self)
% FNAME = WRITEJSIM()

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
[simXMLElem, simXMLDoc] = saveXMLHeader(self, self.model.getLogPath);
[simXMLElem, simXMLDoc] = saveClasses(self, simXMLElem, simXMLDoc);

numOfClasses = length(self.model.classes);
numOfNodes = length(self.model.nodes);
for i=1:(numOfNodes)
    currentNode = self.model.nodes{i,1};
    node = simXMLDoc.createElement('node');
    node.setAttribute('name', currentNode.name);
    
    nodeSections = getSections(currentNode);
    for j=1:length(nodeSections)
        xml_section = simXMLDoc.createElement('section');
        currentSection = nodeSections{1,j};
        if ~isempty(currentSection)
            xml_section.setAttribute('className', currentSection.className);
            switch currentSection.className
                case 'Buffer'
                    xml_section.setAttribute('className', 'Queue'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveBufferCapacity(self, simXMLDoc, xml_section, currentNode);
                    [simXMLDoc, xml_section] = saveDropStrategy(self, simXMLDoc, xml_section);
                    [simXMLDoc, xml_section] = saveGetStrategy(self, simXMLDoc, xml_section, currentNode);
                    [simXMLDoc, xml_section] = savePutStrategy(self, simXMLDoc, xml_section, currentNode);
                case 'Server'
                    [simXMLDoc, xml_section] = saveNumberOfServers(self, simXMLDoc, xml_section, currentNode);
                    [simXMLDoc, xml_section] = saveServerVisits(self, simXMLDoc, xml_section);
                    [simXMLDoc, xml_section] = saveServiceStrategy(self, simXMLDoc, xml_section, currentNode);
                case 'SharedServer'
                    xml_section.setAttribute('className', 'PSServer'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveNumberOfServers(self, simXMLDoc, xml_section, currentNode);
                    [simXMLDoc, xml_section] = saveServerVisits(self, simXMLDoc, xml_section);
                    [simXMLDoc, xml_section] = saveServiceStrategy(self, simXMLDoc, xml_section, currentNode);
                    [simXMLDoc, xml_section] = savePreemptiveStrategy(self, simXMLDoc, xml_section, currentNode);
                    [simXMLDoc, xml_section] = savePreemptiveWeights(self, simXMLDoc, xml_section, currentNode);
                case 'InfiniteServer'
                    xml_section.setAttribute('className', 'Delay'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveServiceStrategy(self, simXMLDoc, xml_section, currentNode);
                case 'LogTunnel'
                    [simXMLDoc, xml_section] = saveLogTunnel(self, simXMLDoc, xml_section, currentNode);
                case 'Dispatcher'
                    xml_section.setAttribute('className', 'Router'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveRoutingStrategy(self, simXMLDoc, xml_section, currentNode);
                case 'StatelessClassSwitcher'
                    xml_section.setAttribute('className', 'ClassSwitch'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveClassSwitchStrategy(self, simXMLDoc, xml_section, currentNode);
                case 'RandomSource'
                    [simXMLDoc, xml_section] = saveArrivalStrategy(self, simXMLDoc, xml_section, currentNode);
                case 'Joiner'
                    xml_section.setAttribute('className', 'Join'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveJoinStrategy(self, simXMLDoc, xml_section, currentNode);
                case 'Forker'
                    xml_section.setAttribute('className', 'Fork'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveForkStrategy(self, simXMLDoc, xml_section, currentNode);
                case 'Storage'
                    xml_section.setAttribute('className', 'Storage'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveTotalCapacity(self, simXMLDoc, xml_section, currentNode);
                    [simXMLDoc, xml_section] = savePlaceCapacities(self, simXMLDoc, xml_section, currentNode);
                    [simXMLDoc, xml_section] = saveDropRule(self, simXMLDoc, xml_section, currentNode);
                    [simXMLDoc, xml_section] = saveGetStrategy(self, simXMLDoc, xml_section, currentNode);
                    [simXMLDoc, xml_section] = savePutStrategies(self, simXMLDoc, xml_section, currentNode);
                case 'Enabling'
                    xml_section.setAttribute('className', 'Enabling'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveEnablingConditions(self, simXMLDoc, xml_section, currentNode);
                    [simXMLDoc, xml_section] = saveInhibitingConditions(self, simXMLDoc, xml_section, currentNode);
                case 'Firing'
                    xml_section.setAttribute('className', 'Firing'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveFiringOutcomes(self, simXMLDoc, xml_section, currentNode);
                case 'Timing'
                    xml_section.setAttribute('className', 'Timing'); %overwrite with JMT class name
                    [simXMLDoc, xml_section] = saveModeNames(self, simXMLDoc, xml_section, currentNode);
                    [simXMLDoc, xml_section] = saveNumbersOfServers(self, simXMLDoc, xml_section, currentNode);
                    [simXMLDoc, xml_section] = saveTimingStrategies(self, simXMLDoc, xml_section, currentNode);
                    [simXMLDoc, xml_section] = saveFiringPriorities(self, simXMLDoc, xml_section, currentNode);
                    [simXMLDoc, xml_section] = saveFiringWeights(self, simXMLDoc, xml_section, currentNode);
            end
            node.appendChild(xml_section);
        end
    end
    simXMLElem.appendChild(node);
end

[simXMLElem, simXMLDoc] = saveMetrics(self, simXMLElem, simXMLDoc);
[simXMLElem, simXMLDoc] = saveLinks(self, simXMLElem, simXMLDoc);

hasReferenceNodes = 0;
preloadNode = simXMLDoc.createElement('preload');
s0 = self.model.getState;
qn = self.model.getStruct;
numOfStations = length(self.model.stations);
for i=1:numOfStations
    isReferenceNode = 0;
    currentNode = self.model.nodes{qn.stationToNode(i),1};
    if (~isa(self.model.stations{i},'Source') && ~isa(self.model.stations{i},'Join'))
        [~, nir] = State.toMarginal(self.model,qn.stationToNode(i),s0{qn.stationToStateful(i)});
        stationPopulationsNode = simXMLDoc.createElement('stationPopulations');
        stationPopulationsNode.setAttribute('stationName', currentNode.name);        
        for r=1:(numOfClasses)
            currentClass = self.model.classes{r,1};
            %        if currentClass.isReferenceStation(currentNode)
            classPopulationNode = simXMLDoc.createElement('classPopulation');
            switch currentClass.type
                case 'open'
                    isReferenceNode = 1;
                    classPopulationNode.setAttribute('population', sprintf('%d',round(nir(r))));
                    classPopulationNode.setAttribute('refClass', currentClass.name);
                    stationPopulationsNode.appendChild(classPopulationNode);                    
                case 'closed'
                    isReferenceNode = 1;
                    %                    classPopulationNode.setAttribute('population', sprintf('%d',currentClass.population));                                        
                    classPopulationNode.setAttribute('population', sprintf('%d',round(nir(r))));
                    classPopulationNode.setAttribute('refClass', currentClass.name);
                    stationPopulationsNode.appendChild(classPopulationNode);
            end
            %        end
        end
    end
    if isReferenceNode
        preloadNode.appendChild(stationPopulationsNode);
    end
    hasReferenceNodes = hasReferenceNodes + isReferenceNode;
end
if hasReferenceNodes
    simXMLElem.appendChild(preloadNode);
end
fname = getJSIMTempPath(self);
try
    xmlwrite(fname, simXMLDoc);
catch
    javaaddpath(which('xercesImpl-2.11.0.jar'));
    javaaddpath(which('xml-apis-2.11.0.jar'));
    pkg load io;
    xmlwrite(fname, simXMLDoc);
end
end