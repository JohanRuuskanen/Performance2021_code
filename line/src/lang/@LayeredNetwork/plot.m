function plot(self,useNodes, showProcs)
% PLOT(SELF,USENODES, SHOWPROCS)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
if ~exist('useNodes','var')
    useNodes = true;
end
if ~exist('showProcs','var')
    showProcs = true;
end

[lqnGraph,taskGraph] = self.getGraph();

%figure;
%plot(taskGraph,'Layout','layered','NodeLabel',taskGraph.Nodes.Node);
%title('Task graph');
figure;
if useNodes
    h = plot(lqnGraph,'Layout','layered','EdgeLabel',lqnGraph.Edges.Weight,'NodeLabel',lqnGraph.Nodes.Node);
else
    h = plot(lqnGraph,'Layout','layered','EdgeLabel',lqnGraph.Edges.Weight,'NodeLabel',lqnGraph.Nodes.Name);
end
title(['Model: ',self.name]);

if showProcs
    for r=findstring(lqnGraph.Nodes.Type,'PS')
        if r>0
            highlight(h,r,'NodeColor','white')
        end
    end
    for r=findstring(lqnGraph.Nodes.Type,'AH')
        if r>0
            highlight(h,r,'NodeColor','cyan')
        end
    end
end

for r=findstring(lqnGraph.Nodes.Type,'T')
    if r>0
        highlight(h,r,'NodeColor','magenta')
    end
end

for r=findstring(lqnGraph.Nodes.Type,'R')
    if r>0
        highlight(h,r,'NodeColor','red')
    end
end

for r=findstring(lqnGraph.Nodes.Type,'H')
    if r>0
        highlight(h,r,'NodeColor','black')
    end
end

for r=findstring(lqnGraph.Nodes.Type,'E')
    if r>0
        highlight(h,r,'NodeColor','green')
    end
end

end
