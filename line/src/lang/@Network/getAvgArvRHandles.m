function [A] = getAvgArvRHandles(self)
% [T] = GETAVGARVRHANDLES()

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

% The method returns the handles to the performance indices but
% they are optional to collect
if isempty(self.handles) || ~isfield(self.handles,'A')
    M = self.getNumberOfStations();
    K = self.getNumberOfClasses();
    
    A = cell(1,K); % arrival rate
    for i=1:M
        for r=1:K
            A{i,r} = Metric(Metric.ArvR, self.classes{r}, self.stations{i});
            %self.addMetric(A{i,r}); % not supported by JMT
            if ~strcmpi(class(self.stations{i}.server),'ServiceTunnel')
                if isempty(self.stations{i}.server.serviceProcess{r}) || strcmpi(class(self.stations{i}.server.serviceProcess{r}{end}),'Disabled')
                    A{i,r}.disable();
                end
            end
        end
    end
    self.handles.A = A;
else
    A = self.handles.A;
end
end
