function task = getNodeTask(self,node)
% TASK = GETNODETASK(SELF,NODE)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
G = self.lqnGraph;
if ischar(node)
    nodeid = findstring(G.Nodes.Name,node);
    task = G.Nodes.Task{nodeid};
else
    task = G.Nodes.Task{node};
end
end
