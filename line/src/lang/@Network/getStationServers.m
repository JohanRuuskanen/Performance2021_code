function S = getStationServers(self)
% S = GETSTATIONSERVERS()

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

for i=1:self.getNumberOfStations()
    S(i,1) = self.stations{i}.numberOfServers;
end
end
