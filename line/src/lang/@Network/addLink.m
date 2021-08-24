function addLink(self, nodeA, nodeB)
% ADDLINK(NODEA, NODEB)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
numberOflinks = length(self.links);
self.links{1+numberOflinks, 1} = {nodeA, nodeB};
end
