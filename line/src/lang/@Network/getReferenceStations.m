function refstat = getReferenceStations(self)
% REFSTAT = GETREFERENCESTATIONS()

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

K = self.getNumberOfClasses();
refstat = zeros(K,1);
for k=1:K
    refstat(k,1) =  findstring(self.getStationNames(),self.classes{k}.reference.name);
end
end
