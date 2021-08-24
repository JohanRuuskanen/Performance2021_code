function entry = getNodeEntry(self,node)
% ENTRY = GETNODEENTRY(SELF,NODE)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
G = self.lqnGraph;
if ischar(node)
    nodeid = findstring(G.Nodes.Name,node);
    entry = G.Nodes.Entry{nodeid};
else
    entry = G.Nodes.Entry{node};
end
end
