function [sched, schedid, schedparam] = refreshScheduling(self)
% [SCHED, SCHEDID, SCHEDPARAM] = REFRESHSCHEDULING()
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

% determine scheduling parameters
M = self.getNumberOfStations();
K = self.getNumberOfClasses();

sched = self.getStationScheduling();
schedparam = zeros(M,K);
for i=1:M
    if isempty(self.getIndexSourceStation) || i ~= self.getIndexSourceStation
        switch self.stations{i}.server.className
            case 'ServiceTunnel'
                % do nothing
            otherwise
                if ~isempty(self.stations{i}.schedStrategyPar) & ~isnan(self.stations{i}.schedStrategyPar) %#ok<AND2>
                    schedparam(i,:) = self.stations{i}.schedStrategyPar;
                else
                    switch sched(i)
                        case SchedStrategy.SEPT
                            svcTime = zeros(1,K);
                            for k=1:K
                                svcTime(k) = self.nodes{i}.serviceProcess{k}.getMean;
                            end
                            [svcTimeSorted] = sort(unique(svcTime));
                            self.nodes{i}.schedStrategyPar = zeros(1,K);
                            for k=1:K
                                self.nodes{i}.schedStrategyPar(k) = find(svcTimeSorted == svcTime(k));
                            end
                        case SchedStrategy.LEPT
                            svcTime = zeros(1,K);
                            for k=1:K
                                svcTime(k) = self.nodes{i}.serviceProcess{k}.getMean;
                            end
                            [svcTimeSorted] = sort(unique(svcTime),'descend');
                            self.nodes{i}.schedStrategyPar = zeros(1,K);
                            for k=1:K
                                self.nodes{i}.schedStrategyPar(k) = find(svcTimeSorted == svcTime(k));
                            end
                    end
                end
        end
    end
end

if ~isempty(self.qn)
    self.qn.setSched(sched, schedparam);
end
end
