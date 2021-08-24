function self = initDefault(self)
% SELF = INITDEFAULT()

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
self.lqnGraph.Edges.Weight=log(self.lqnGraph.Edges.Weight); % use log-weights
V=self.lqnGraph.distances; V=exp(V); V(isinf(V))=0; % mean number of calls to the node from ref task
D=self.lqnGraph.Nodes.D;
D(~isfinite(D))=0;
D(isnan(D))=0;
respTime = V*D; % this is for entries activities with calls
self.param.Nodes.RespT = respTime; % response time for one visit
self.param.Nodes.Util = (self.lqnGraph.Nodes.D)./respTime;
self.param.Nodes.Tput = 1./respTime;
QLen = self.param.Nodes.Tput .* self.param.Nodes.RespT;
self.param.Nodes.Util(findstring(self.lqnGraph.Nodes.Type,'E')) = QLen(findstring(self.lqnGraph.Nodes.Type,'E'));
self.lqnGraph.Edges.Weight=exp(self.lqnGraph.Edges.Weight); % restore initial weights
self.param.Edges.RespT(:)=0; % restore initial weights
for e=find(self.lqnGraph.Edges.Type==1)' % initialize resp time of sync calls to demand values
    self.param.Edges.RespT(e) = self.param.Nodes.RespT(self.getNodeIndex(self.lqnGraph.Edges.EndNodes{e,2}));
    self.param.Edges.Tput(e) = self.param.Nodes.Tput(self.getNodeIndex(self.lqnGraph.Edges.EndNodes{e,2}));
end
end
