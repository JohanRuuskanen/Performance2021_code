function addNode(self, node)
% ADDNODE(NODE)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

%
if sum(cellfun(@(x) strcmp(x.name,node.name), {self.nodes{1:end}}))>0
    line_error(mfilename,'A node with an identical name already exists.');
end
self.nodes{end+1,1} = node;
if isa(node,'Station')
    self.stations{end+1,1} = node;
end
self.setUsedFeatures(class(node)); % station type
self.setUsedFeatures(class(node.input)); % station type
self.setUsedFeatures(class(node.server)); % station type
self.setUsedFeatures(class(node.output)); % station type
end
