function refreshStruct(self)
% REFRESHSTRUCT()
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
self.sanitize();
nodetypes = self.getNodeTypes();
classnames = self.getClassNames();
nodenames = self.getNodeNames();
jobs = self.getNumberOfJobs();
servers = self.getStationServers();
refstat = self.getReferenceStations();
routing = self.getRoutingStrategies();
self.qn = NetworkStruct(nodetypes, nodenames, classnames, servers, jobs(:), refstat, routing);
self.refreshPriorities();
self.refreshService();
wantVisits = true;
if any(nodetypes == NodeType.Cache)
    wantVisits = false;
end
self.refreshChains(wantVisits);
self.refreshLocalVars(); % depends on chains (rtnodes)
self.refreshSync(); % this assumes that refreshChain is called before
%self.qn.forks = self.getForks(self.qn.rt);
end
