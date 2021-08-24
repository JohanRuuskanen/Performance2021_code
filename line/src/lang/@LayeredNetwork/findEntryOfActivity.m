function [entryName, entryFullName] = findEntryOfActivity(self,activity)
% [ENTRYNAME, ENTRYFULLNAME] = FINDENTRYOFACTIVITY(SELF,ACTIVITY)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

G = self.lqnGraph;
actIdx = self.getNodeIndex(activity);
% we need to search the graph
prec = G.predecessors(actIdx);
if isempty(prec)
    entryName = '';
    entryFullName = '';
    return
end
p = prec(1);
if strcmpi(G.Nodes.Type{p},'A')
    % if it is an activity, go recursively
    [entryName, entryFullName] = self.findEntryOfActivity(G.Nodes.Name{p});
else
    % otherwise, we found the entry
    entryName = G.Nodes.Name{p};
    entryFullName = G.Nodes.Node{p};
end
end
