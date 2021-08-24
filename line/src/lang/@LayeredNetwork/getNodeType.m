function type = getNodeType(self,nodeid)
% TYPE = GETNODETYPE(SELF,NODEID)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
G = self.lqnGraph;
if ischar(nodeid)
    type = G.Nodes.Type{self.getNodeIndex(nodeid)};
else
    type = G.Nodes.Type{nodeid};
end
end
