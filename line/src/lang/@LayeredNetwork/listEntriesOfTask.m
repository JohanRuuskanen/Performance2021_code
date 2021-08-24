function entries=listEntriesOfTask(self,task)
% ENTRIES=LISTENTRIESOFTASK(SELF,TASK)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
G = self.lqnGraph;
entries = {};
if ischar(task)
    taskid = self.getNodeIndex(task);
else
    taskid = task;
end
for s=G.successors(taskid)'
    name=self.getNodeName(s);
    if name(1)=='E'
        entries{end+1}=name;
    end
end
end

