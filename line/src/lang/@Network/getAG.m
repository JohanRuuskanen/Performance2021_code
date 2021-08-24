% getAG : export model in agent representation
function ag = getAG(self)
% AG = GETAG()

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

% parses all but the service processes
if isempty(self.ag)
    self.refreshAG();
end
ag=self.ag;
end
