function idx = findEdgeIndex(self,source,dest)
% IDX = FINDEDGEINDEX(SELF,SOURCE,DEST)
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

%G = self.lqnGraph;
%idx = G.findedge(source, dest);
%return
if ischar(source)
    idx = matchrow(self.endNodes,[self.getNodeIndex(source),self.getNodeIndex(dest)]);
else
    idx = matchrow(self.endNodes,[source,dest]);
end
%% code below is correct but slow
% if ischar(source)
%     from = findstring(G.Edges.EndNodes(:,1),source);
%     to = findstring(G.Edges.EndNodes(:,2),dest);
% else
%     from = findstring(G.Edges.EndNodes(:,1),self.getNodeName(source));
%     to = findstring(G.Edges.EndNodes(:,2),self.getNodeName(dest));
% end
% idx = intersect(from,to);
end
