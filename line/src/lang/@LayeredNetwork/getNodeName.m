function name = getNodeName(self,node,useNode)
% NAME = GETNODENAME(SELF,NODE,USENODE)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
G = self.lqnGraph;
if ~exist('useNode','var') || useNode == false
    name = G.Nodes.Name{node};
else
    name = G.Nodes.Node{node};
end
end
