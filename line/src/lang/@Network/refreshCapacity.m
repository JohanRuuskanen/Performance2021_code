function [capacity, classcap] = refreshCapacity(self)
% [CAPACITY, CLASSCAP] = REFRESHCAPACITY()

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
M = self.getNumberOfStations();
K = self.getNumberOfClasses();
C = self.qn.nchains;
% set zero buffers for classes that are disabled
classcap = Inf*ones(M,K);
chaincap = Inf*ones(M,C);
capacity = zeros(M,1);
for i=1:M
    station = self.getStationByIndex(i);
    for r=1:K
        if isa(station, 'Place')
            classcap(i,r) = min(station.classCap(r), station.cap);
        elseif isempty(self.qn.rates(i,r)) || self.qn.rates(i,r)==0 || any(~isfinite(self.qn.rates(i,r)))
            classcap(i,r) = 0;
            chaincap(i,r) = 0;
        else
            c = find(self.qn.chains(:,r),1,'first'); % chain of class r
            chaincap(i,c) = sum(self.qn.njobs(find(self.qn.chains(c,:)))); %#ok<FNDSB>
            classcap(i,r) = sum(self.qn.njobs(find(self.qn.chains(c,:)))); %#ok<FNDSB>
            if station.classCap(r) >= 0
                classcap(i,r) = min(classcap(i,r), station.classCap(r));
            end
            if station.cap >= 0
                classcap(i,r) = min(classcap(i,r),station.cap);
            end
        end
    end
    capacity(i,1) = sum(chaincap(i,:));
end
if ~isempty(self.qn) %&& isprop(self.qn,'cap')
    self.qn.setCapacity(capacity, classcap);
end
end
