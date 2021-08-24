function submodels = getGraphLayers(self)
% SUBMODELS = GETGRAPHLAYERS()

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

submodels = layerize_loose(self);
end

function submodels = layerize_loose(self)
% SUBMODELS = LAYERIZE_LOOSE(SELF)

% first build all layers
reftasks = findstring(self.taskGraph.Nodes.Type,'R');
order = toposort(self,reftasks);
order = order(~ismember(order,reftasks));
submodels = cell(1,length(order));
for l = 1:length(submodels)
    submodels{l} = digraph();
    nameTo = self.taskGraph.Nodes.Name{order(l)};
    submodels{l} = submodels{l}.addnode(nameTo);
    pred = predecessors(self.taskGraph,order(l));
    for p = pred'
        nameFrom = self.taskGraph.Nodes.Name{p};
        submodels{l} = submodels{l}.addedge(nameFrom,nameTo);
    end
end

% update full name
for l = 1:length(submodels)
    submodels{l}.Nodes.Node = cell(submodels{l}.numnodes,1);
    for i = 1:submodels{l}.numnodes
        j = self.lqnGraph.findnode(submodels{l}.Nodes.Name(i));
        submodels{l}.Nodes.Node(i) = self.lqnGraph.Nodes.Node(j);
    end
end
end

function order = toposort(self,reftasks)
% ORDER = TOPOSORT(SELF,REFTASKS)

self.topoOrder = [];
self.dfsMarks = zeros(height(self.taskGraph.Nodes),1);
for r = reftasks'
    dfsearch(self,r);
end
self.topoOrder = flip(self.topoOrder);
order = self.topoOrder;
end

function dfsearch(self,n)
% DFSEARCH(SELF,N)

if self.dfsMarks(n) > 0
    return
end
self.dfsMarks(n) = 1;
succ = successors(self.taskGraph,n);
for s = succ'
    dfsearch(self,s);
end
self.dfsMarks(n) = 2;
self.topoOrder(end+1,1) = n;
end
