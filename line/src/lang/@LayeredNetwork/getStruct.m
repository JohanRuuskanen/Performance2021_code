function lqn = getStruct(self)
% LQN = GETSTRUCT(SELF)
%
%
% Copyright 2012-2020, Imperial College London

lqn = struct();
lqn.nidx = 0;
lqn.nhosts = length(self.hosts);
lqn.ntasks = length(self.tasks);
lqn.nreftasks = length(self.reftasks);
lqn.nacts = length(self.activities);
lqn.nentries = length(self.entries);
lqn.ntasksof = [];
lqn.nentriesof = [];
lqn.nactsof = [];
lqn.tshift = lqn.nhosts;
lqn.eshift = lqn.nhosts + lqn.ntasks;
lqn.ashift = lqn.nhosts + lqn.ntasks + lqn.nentries;

for p=1:lqn.nhosts
    lqn.ntasksof(p) = length(self.hosts{p}.tasks);
end

%% analyze static properties
lqn.nidx = lqn.nhosts + lqn.ntasks + lqn.nentries + lqn.nacts;
idx = 1;
lqn.hostidx = [];
lqn.taskidx = [];
lqn.entryidx = [];
lqn.actidx = [];
lqn.tasksof = cell(lqn.nhosts,1);
lqn.entriesof = cell(lqn.ntasks,1);
lqn.actsof = cell(lqn.ntasks,1);
lqn.callsof = cell(lqn.ntasks,1);
lqn.hostdem = {};
lqn.think = {};
lqn.sched = categorical([]);
lqn.names = {};
lqn.hashnames = {};
lqn.shortnames = {};
lqn.mult = zeros(lqn.nhosts+lqn.ntasks,1);
lqn.repl = zeros(lqn.nhosts+lqn.ntasks,1);
lqn.type = zeros(lqn.nidx,1);
lqn.graph = sparse([]);
lqn.replies = [];
lqn.replygraph = sparse([]);

tshift = lqn.nhosts;
eshift = lqn.nhosts + lqn.ntasks;
ashift = lqn.nhosts + lqn.ntasks + lqn.nentries;

lqn.parent = [];
for p=1:lqn.nhosts
    lqn.hostidx(end+1) = idx;
    lqn.sched(idx,1) = SchedStrategy.fromText(self.hosts{p}.scheduling);
    lqn.mult(idx,1) = self.hosts{p}.multiplicity;
    lqn.repl(idx,1) = self.hosts{p}.replication;
    lqn.names{idx,1} = self.hosts{p}.name;
    lqn.hashnames{idx,1} = ['P:',lqn.names{idx,1}];
    lqn.shortnames{idx,1} = ['P',num2str(p)];
    lqn.type(idx,1) = LayeredNetworkElement.HOST; % processor
    idx = idx + 1;
end

for p=1:lqn.nhosts
    pidx = p;
    for t=1:lqn.ntasksof(p)
        lqn.taskidx(end+1) = idx;
        lqn.sched(idx,1) = SchedStrategy.fromText(self.hosts{p}.tasks(t).scheduling);
        lqn.hostdem{idx,1} = Immediate();
        lqn.think{idx,1} = self.hosts{p}.tasks(t).thinkTime;
        lqn.mult(idx,1) = self.hosts{p}.tasks(t).multiplicity;
        lqn.repl(idx,1) = self.hosts{p}.tasks(t).replication;
        lqn.names{idx,1} = self.hosts{p}.tasks(t).name;        
        switch lqn.sched(idx,1)
            case SchedStrategy.REF
                lqn.hashnames{idx,1} = ['R:',lqn.names{idx,1}];
                lqn.shortnames{idx,1} = ['R',num2str(idx-tshift)];
            otherwise
                lqn.hashnames{idx,1} = ['T:',lqn.names{idx,1}];
                lqn.shortnames{idx,1} = ['T',num2str(idx-tshift)];
        end
        lqn.parent(idx) = pidx;
        lqn.graph(idx, pidx) = 1;
        lqn.nentriesof(idx) = length(self.hosts{p}.tasks(t).entries);
        lqn.nactsof(idx) = length(self.hosts{p}.tasks(t).activities);
        lqn.type(idx) = LayeredNetworkElement.TASK; % task
        idx = idx + 1;
    end
    lqn.tasksof{pidx} = find(lqn.parent == pidx);
end


tasks = self.tasks;
for t = 1:lqn.ntasks
    tidx = lqn.taskidx(t);
    for e=1:lqn.nentriesof(tidx)
        lqn.entryidx(end+1) = idx;
        lqn.names{idx,1} = self.tasks{t}.entries(e).name;
        lqn.hashnames{idx,1} = ['E:',lqn.names{idx,1}];
        lqn.shortnames{idx,1} = ['E',num2str(idx-eshift)];
        lqn.hostdem{idx,1} = Immediate();
        lqn.parent(idx) = tidx;
        lqn.graph(tidx,idx) = 1;
        lqn.entriesof{tidx}(e) = idx;
        lqn.type(idx) = LayeredNetworkElement.ENTRY; % entries
        idx = idx + 1;
    end
end

for t = 1:lqn.ntasks
    tidx = lqn.taskidx(t);
    for a=1:lqn.nactsof(tidx)
        lqn.actidx(end+1) = idx;
        lqn.names{idx,1} = tasks{t}.activities(a).name;
        lqn.hashnames{idx,1} = ['A:',lqn.names{idx,1}];
        lqn.shortnames{idx,1} = ['A',num2str(idx - ashift)];
        lqn.hostdem{idx,1} = tasks{t}.activities(a).hostDemand;
        lqn.parent(idx) = tidx;
        lqn.actsof{tidx}(a) = idx;
        lqn.type(idx) = LayeredNetworkElement.ACTIVITY; % activities
        idx = idx + 1;
    end
end

nidx = idx - 1; % number of indices
lqn.graph(nidx,nidx) = 0;

%% now analyze calls
cidx = 0;
lqn.callidx = sparse(lqn.nidx,lqn.nidx);
lqn.calltype = sparse([],lqn.nidx,1);
lqn.iscaller = sparse(lqn.nidx,lqn.nidx);
lqn.issynccaller = sparse(lqn.nidx,lqn.nidx);
lqn.isasynccaller = sparse(lqn.nidx,lqn.nidx);
lqn.callpair = [];
lqn.callproc = {};
lqn.callnames = {};
lqn.callhashnames = {};
lqn.callshortnames = {};
lqn.taskgraph = sparse([],lqn.ntasks, lqn.ntasks);
lqn.actpre = sparse(lqn.nidx,1);
lqn.actpost = sparse(lqn.nidx,1);

for t = 1:lqn.ntasks
    tidx = lqn.taskidx(t);
    lqn.actsof{tidx} = zeros(1,lqn.nactsof(tidx));
    for a=1:lqn.nactsof(tidx)
        aidx = findstring(lqn.hashnames, ['A:',tasks{t}.activities(a).name]);
        lqn.callsof{aidx} = [];
        lqn.actsof{tidx}(a) = aidx;        
        boundToEntry = tasks{t}.activities(a).boundToEntry;
        %for b=1:length(boundToEntry)
        eidx = findstring(lqn.hashnames, ['E:',boundToEntry]);
        if eidx>0
            lqn.graph(eidx, aidx) = 1;
        end
        %end
        
        for s=1:length(tasks{t}.activities(a).syncCallDests)
            target_eidx = findstring(lqn.hashnames, ['E:',tasks{t}.activities(a).syncCallDests{s}]);
            target_tidx = lqn.parent(target_eidx);
            cidx = cidx + 1;
            lqn.callidx(aidx, target_eidx) = cidx;
            lqn.calltype(cidx,1) = CallType.ID_SYNC;
            lqn.callpair(cidx,1:2) = [aidx,target_eidx];
            lqn.callnames{cidx,1} = [lqn.names{aidx},'=>',lqn.names{target_eidx}];
            lqn.callhashnames{cidx,1} = [lqn.hashnames{aidx},'=>',lqn.hashnames{target_eidx}];
            lqn.callshortnames{cidx,1} = [lqn.shortnames{aidx},'=>',lqn.shortnames{target_eidx}];
            lqn.callproc{cidx,1} = Geometric(1/tasks{t}.activities(a).syncCallMeans(s)); % synch
            lqn.callsof{aidx}(end+1) = cidx;
            lqn.iscaller(tidx, target_tidx) = true;
            lqn.iscaller(tidx, target_eidx) = true;
            lqn.iscaller(aidx, target_eidx) = true;
            lqn.issynccaller(tidx, target_tidx) = true;
            lqn.issynccaller(tidx, target_eidx) = true;
            lqn.issynccaller(aidx, target_eidx) = true;
            lqn.taskgraph(tidx, target_tidx) = 1;
            lqn.graph(aidx, target_eidx) = 1;
        end
        
        for s=1:length(tasks{t}.activities(a).asyncCallDests)
            target_eidx = findstring(lqn.hashnames,['E:',tasks{t}.activities(a).asyncCallDests{s}]);
            target_tidx = lqn.parent(target_eidx);
            cidx = cidx + 1;
            lqn.callidx(aidx, target_eidx) = cidx;
            lqn.callpair(cidx,1:2) = [aidx,target_eidx];
            lqn.calltype(cidx,1) = CallType.ID_ASYNC; % async
            lqn.callnames{cidx,1} = [lqn.names{aidx},'->',lqn.names{target_eidx}];
            lqn.callhashnames{cidx,1} = [lqn.hashnames{aidx},'->',lqn.hashnames{target_eidx}];
            lqn.callshortnames{cidx,1} = [lqn.shortnames{aidx},'->',lqn.shortnames{target_eidx}];
            lqn.callproc{cidx,1} = Geometric(1/tasks{t}.activities(a).asyncCallMeans(s)); % asynch
            lqn.callsof{aidx}(end+1) = cidx;
            lqn.iscaller(tidx, target_tidx) = true;
            lqn.iscaller(tidx, target_eidx) = true;
            lqn.isasynccaller(tidx, target_tidx) = true;
            lqn.isasynccaller(tidx, target_eidx) = true;
            lqn.isasynccaller(aidx, target_eidx) = true;
            lqn.taskgraph(tidx, target_tidx) = 1;
            lqn.graph(aidx, target_eidx) = 1;
        end
    end
    
    for ap=1:length(tasks{t}.precedences)
        %pretype = tasks{t}.precedences(ap).preType;
        preacts = tasks{t}.precedences(ap).preActs;
        postacts = tasks{t}.precedences(ap).postActs;
        for prea = 1:length(preacts)
            preaidx = findstring(lqn.hashnames, ['A:',tasks{t}.precedences(ap).preActs{prea}]);
            if isempty(tasks{t}.precedences(ap).preParams)
                preParam = 1.0; 
            else
                preParam = tasks{t}.precedences(ap).preParams(prea); 
            end
            for posta = 1:length(postacts)
                postaidx = findstring(lqn.hashnames, ['A:',tasks{t}.precedences(ap).postActs{posta}]);                
                if isempty(tasks{t}.precedences(ap).postParams)
                    postParam = 1.0; 
                else
                    postParam = tasks{t}.precedences(ap).postParams(posta); 
                end
                lqn.graph(preaidx, postaidx) = preParam * postParam;
                lqn.actpre(preaidx) = sparse(ActivityPrecedence.getPrecedenceId(tasks{t}.precedences(ap).preType));
                lqn.actpost(preaidx) = sparse(ActivityPrecedence.getPrecedenceId(tasks{t}.precedences(ap).postType));
            end
        end
    end
end

lqn.replies = false(1,lqn.nidx);
lqn.replygraph = 0*lqn.graph;
for t = 1:lqn.ntasks
    tidx = lqn.taskidx(t);
    for aidx = lqn.actsof{tidx}
        postaidxs = find(lqn.graph(aidx, :));
        isreply = true;
        % if no successor is an action of tidx
        for postaidx = postaidxs
            if any(lqn.actsof{tidx} == postaidx)
                isreply = false;
            end
        end
        if isreply
            % this is a leaf node, search backward for the parent entry,
            % which is assumed to be unique
            lqn.replies(aidx) = true;
            parentidx = aidx;
            while lqn.type(parentidx) ~= LayeredNetworkElement.ENTRY
                ancestors = find(lqn.graph(:,parentidx));
                parentidx = at(ancestors,1); % only choose first ancestor
            end
            if lqn.type(parentidx) == LayeredNetworkElement.ENTRY
                lqn.replygraph(aidx, parentidx) = 1;
            end
        end
    end
end
lqn.ncalls = size(lqn.calltype,1);

% correct multiplicity for infinite server stations
for tidx = find(lqn.sched== SchedStrategy.INF)
    if lqn.type(tidx) == LayeredNetworkElement.TASK
        callers = find(lqn.taskgraph(:, tidx));
        callers_inf = strcmp(lqn.mult(callers), SchedStrategy.INF);
        if any(callers_inf) 
            % if a caller is also inf, then we would need to recursively
            % determine the maximum multiplicity, we instead use a
            % heuristic value
            lqn.mult(tidx) = sum(lqn.mult(~callers_inf)) + sum(callers_inf)*max(lqn.mult);
        else
            lqn.mult(tidx) = sum(lqn.mult(callers));
        end
    end
end

lqn.isref = lqn.sched==SchedStrategy.REF;
end