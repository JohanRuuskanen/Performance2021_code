function idx = getNodeIndexInTaskGraph(self,node)
% IDX = GETNODEINDEXINTASKGRAPH(SELF,NODE)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
G = self.taskGraph;
idx = findstring(G.Nodes.Name,node);
%idx = H.findnode(node);
end
