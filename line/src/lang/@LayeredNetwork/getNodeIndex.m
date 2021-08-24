function idx = getNodeIndex(self,node)
% IDX = GETNODEINDEX(SELF,NODE)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

%G = self.lqnGraph;
%idx = findstring(G.Nodes.Name,node);
%idx = G.findnode(node);
idx = findstring(self.nodeNames,node);
end
