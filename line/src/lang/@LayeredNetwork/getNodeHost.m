function proc = getNodeHost(self,node)
% PROC = GETNODEHOST(SELF,NODE)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
%G = self.lqnGraph;
if ischar(node)
    id = findstring(self.nodeNames,node);
    %proc = G.Nodes.Proc{id};
    proc = self.nodeNames{self.nodeDep(id,1)};
else % index
    %proc = G.Nodes.Proc{node};
    proc = self.nodeNames{self.nodeDep(node,1)};
end
end
