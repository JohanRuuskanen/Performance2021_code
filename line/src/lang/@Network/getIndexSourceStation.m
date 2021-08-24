function sourceidx = getIndexSourceStation(self)
% INDEX = GETINDEXSOURCESTATION()

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
if isempty(self.sourceidx)
    self.sourceidx = find(cellisa(self.stations,'Source'));
end
sourceidx = self.sourceidx;
end
