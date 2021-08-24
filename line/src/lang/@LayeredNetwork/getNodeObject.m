function obj = getNodeObject(self,node)
% OBJ = GETNODEOBJECT(SELF,NODE)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
G = self.lqnGraph;
if ischar(node)
    obj = G.Nodes.Object{self.getNodeIndex(node)};
else
    obj = G.Nodes.Object{node};
end
end
