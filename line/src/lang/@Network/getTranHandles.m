%% getTranHandles: add all transient mean performance indexes
% Q{i,r}: timeseries of mean queue-length of class r at node i
% U{i,r}: timeseries of mean utilization of class r at node i
% R{i,r}: timeseries of mean response time of class r at node i (summed across visits)
% T{i,r}: timeseries of mean throughput of class r at node i
function [Qt,Ut,Tt] = getTranHandles(self)
% [QT,UT,TT] = GETTRANHANDLES()

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

% The method returns the handles to the performance indices but
% they are optional to collect
M = self.getNumberOfStations();
K = self.getNumberOfClasses();

Tt = cell(1,K); % throughputs
Qt = cell(M,K); % queue-length
%Rt = cell(M,K); % response times
Ut = cell(M,1); % utilizations
for i=1:M
    for r=1:K
        Tt{i,r} = Metric(Metric.TranTput, self.classes{r}, self.stations{i});
        self.addTranMetric(Tt{i,r});
        Qt{i,r} = Metric(Metric.TranQLen, self.classes{r}, self.stations{i});
        self.addTranMetric(Qt{i,r});
        %        Rt{i,r} = Metric(Metric.TranRespT, self.classes{r}, self.stations{i});
        %        self.addTranMetric(Rt{i,r});
        Ut{i,r} = Metric(Metric.TranUtil, self.classes{r}, self.stations{i});
        self.addTranMetric(Ut{i,r});
        if isa(self.stations{i},'Source')
            Qt{i,r}.disable();
            %            Rt{i,r}.disable();
            Ut{i,r}.disable();
        end
        if isa(self.stations{i},'Sink')
            Qt{i,r}.disable();
            %            Rt{i,r}.disable();
            Ut{i,r}.disable();
        end
        if isa(self.stations{i},'Join') || isa(self.stations{i},'Fork')
            Ut{i,r}.disable();
        end
        if ~strcmpi(class(self.stations{i}.server),'ServiceTunnel')
            if isempty(self.stations{i}.server.serviceProcess{r}) || strcmpi(class(self.stations{i}.server.serviceProcess{r}{end}),'Disabled')
                Tt{i,r}.disable();
                Qt{i,r}.disable();
                %                Rt{i,r}.disable();
                Ut{i,r}.disable();
            end
        end
    end
end
end
