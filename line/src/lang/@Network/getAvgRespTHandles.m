% R(i,r): mean response time of class r at node i (summed across visits)
function R = getAvgRespTHandles(self)
% R = GETAVGRESPTHANDLES()

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

% The method returns the handles to the performance indices but
% they are optional to collect
if isempty(self.handles) || ~isfield(self.handles,'R')
    M = self.getNumberOfStations();
    K = self.getNumberOfClasses();
    
    R = cell(M,K); % response times
    for i=1:M
        for r=1:K
            R{i,r} = Metric(Metric.RespT, self.classes{r}, self.stations{i});
            if isa(self.stations{i},'Source')
                R{i,r}.disable();
            end
            if isa(self.stations{i},'Sink')
                R{i,r}.disable();
            end
            if ~strcmpi(class(self.stations{i}.server),'ServiceTunnel')
                if isempty(self.stations{i}.server.serviceProcess{r}) || strcmpi(class(self.stations{i}.server.serviceProcess{r}{end}),'Disabled')
                    R{i,r}.disable();
                end
            end
        end
    end
    self.addMetric(R);
    self.handles.R = R;
else
    R = self.handles.R;
end
end
