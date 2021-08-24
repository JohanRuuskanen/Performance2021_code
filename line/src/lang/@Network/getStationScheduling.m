function sched = getStationScheduling(self)
% SCHED = GETSTATIONSCHEDULING()

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

for i=1:self.getNumberOfStations()
    if isinf(self.stations{i}.numberOfServers)
        sched(i,1) = SchedStrategy.INF;
    else
        if i == self.getIndexSourceStation()
            sched(i,1) = SchedStrategy.EXT;
        else
            sched(i,1) = self.stations{i}.schedStrategy;
        end
    end
end
end
