function [simElem,simDoc] = saveXMLHeader(self, logPath)
% [SIMELEM,SIMDOC] = SAVEXMLHEADER(LOGPATH)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
fname = [self.getFileName(), ['.', 'jsimg']];
if isoctave
    try
        simDoc = javaObject('org.apache.xerces.dom.DocumentImpl');
        simElem = simDoc.createElement('sim');
    catch
        javaaddpath(which('xercesImpl-2.11.0.jar'));
        javaaddpath(which('xml-apis-2.11.0.jar'));
        pkg load io;
        simDoc = javaObject('org.apache.xerces.dom.DocumentImpl');
        simElem = simDoc.createElement('sim');
        simDoc.appendChild(simElem);
    end
else
    simDoc = com.mathworks.xml.XMLUtils.createDocument('sim');
    simElem = simDoc.getDocumentElement;
end
simElem.setAttribute('xmlns:xsi', self.xmlnsXsi);
simElem.setAttribute('name', fname);
%simElem.setAttribute('timestamp', '"Tue Jan 1 00:00:01 GMT+00:00 2000"');
simElem.setAttribute('xsi:noNamespaceSchemaLocation', 'SIMmodeldefinition.xsd');
simElem.setAttribute('disableStatisticStop', 'true');
simElem.setAttribute('logDecimalSeparator', '.');
simElem.setAttribute('logDelimiter', ';');
simElem.setAttribute('logPath', logPath);
simElem.setAttribute('logReplaceMode', '0');
simElem.setAttribute('maxSamples', int2str(self.maxSamples));
simElem.setAttribute('maxEvents', int2str(self.maxEvents));
if ~isinf(self.maxSimulatedTime)
    simElem.setAttribute('maxSimulated', num2str(self.maxSimulatedTime,'%.3f'));
end
simElem.setAttribute('polling', '1.0');
simElem.setAttribute('seed', int2str(self.options.seed));
end
