function [self]=setGraph(self,lqnGraph,taskGraph)
% [SELF]=SETGRAPH(SELF,LQNGRAPH,TASKGRAPH)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
self.lqnGraph = lqnGraph;
if nargin>2
    self.taskGraph = taskGraph;
end
end

