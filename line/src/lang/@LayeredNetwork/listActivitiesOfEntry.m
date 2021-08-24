function acts=listActivitiesOfEntry(self,entry)
% ACTS=LISTACTIVITIESOFENTRY(SELF,ENTRY)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
G = self.lqnGraph;
acts = G.Nodes.Name(findstring(G.Nodes.Entry,entry));
acts = {acts{:}};
end

