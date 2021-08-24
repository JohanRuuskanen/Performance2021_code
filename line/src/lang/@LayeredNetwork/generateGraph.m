function self = generateGraph(self)
% SELF = GENERATEGRAPH()

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

self.lqnGraph = digraph();
% add hosts
self.nodeDep = [];
fullname = {};
name = {};
type = {};
object = {};
EndNodes = [];
Weight = [];
Object = [];
demand = [];
multiplicity = [];
%boundToEntry = [];
host_proc = {};
host_task = {};
host_entry = {};
maxjobs = [];
nodeDep = [];
ctrp = 0; ctrt = 0; ctre = 0; ctrap = 0; ctras = 0;
for p=1:length(self.hosts)
    ctrp = ctrp + 1;
    fullname{end+1} = self.hosts{p}.name;
    type{end+1} = 'H';
    hostname = sprintf('P%d',ctrp);
    name{end+1} = sprintf('P%d',ctrp);
    pidx = length(name);
    object{end+1} = self.hosts{p};
    demand(end+1) = 0.0;
    multiplicity(end+1)= object{end}.multiplicity;
    %    boundToEntry(end+1)=false;
    host_proc{end+1}= '';
    host_task{end+1}= '';
    host_entry{end+1}= '';
    nodeDep(end+1,1:3) = [NaN,NaN,NaN];
    for t=1:length(self.hosts{p}.tasks)
        ctrt = ctrt + 1;
        fullname{end+1} = self.hosts{p}.tasks(t).name;
        name{end+1} = sprintf('T%d',ctrt);
        tidx = length(name);
        taskname = name{end};
        demand(end+1) = self.hosts{p}.tasks(t).thinkTimeMean;
        object{end+1} = self.hosts{p}.tasks(t);
        multiplicity(end+1)= object{end}.multiplicity;
        %        boundToEntry(end+1)=false;
        host_proc{end+1}= hostname;
        host_task{end+1}= '';
        host_entry{end+1}= '';
        nodeDep(end+1,1:3) = [pidx,NaN,NaN];
        switch self.hosts{p}.tasks(t).scheduling
            case 'ref'
                type{end+1} = 'R'; % reference task
                maxjobs(:,end+1) =  multiplicity(end);
            otherwise
                type{end+1} = 'T';
        end
        for e=1:length(self.hosts{p}.tasks(t).entries)
            ctre = ctre + 1;
            fullname{end+1} = self.hosts{p}.tasks(t).entries(e).name;
            name{end+1} = sprintf('E%d',ctre);
            eidx = length(name);
            entryname = name{end};
            type{end+1} = 'E';
            demand(end+1) = 0.0;
            object{end+1} = self.hosts{p}.tasks(t).entries(e);
            multiplicity(end+1)= NaN;
            %            boundToEntry(end+1)=false;
            host_proc{end+1}= hostname;
            host_task{end+1}= taskname;
            host_entry{end+1}= '';
            nodeDep(end+1,1:3) = [pidx,tidx,NaN];
        end
        for a=1:length(self.hosts{p}.tasks(t).activities)
            ctras = ctras + 1;
            fullname{end+1} = self.hosts{p}.tasks(t).activities(a).name;
            demand(end+1) = self.hosts{p}.tasks(t).activities(a).hostDemandMean;
            %            name{end+1} = sprintf('A%d(%.3f)',ctras,demand(end));
            name{end+1} = sprintf('A%d',ctras);
            type{end+1} = 'A';
            object{end+1} = self.hosts{p}.tasks(t).activities(a);
            multiplicity(end+1)= NaN;
            %            boundToEntry(end+1)=false;
            %            if ~isempty(object{end}.boundToEntry)
            %                boundToEntry(end)=true;
            %            end
            host_proc{end+1}= hostname;
            host_task{end+1}= taskname;
            host_entry{end+1}= 'NaN';
            if ~isempty(object{end}.boundToEntry)
                host_entry{end}= name{findstring(fullname,object{end}.boundToEntry)};
                eidx = findstring(name,host_entry{end});
                nodeDep(end+1,1:3) = [pidx,tidx,eidx];
            else
                nodeDep(end+1,1:3) = [pidx,tidx,NaN];
            end
        end
    end
end

myTable = Table();
myTable.Name = name(:);
self.nodeNames = name(:);
myTable.Type = type(:);
myTable.Proc = host_proc(:);
myTable.Task = host_task(:);
myTable.Entry = host_entry(:);
myTable.D = demand(:);
myTable.Mult = multiplicity(:);
%myTable.BtoE = boundToEntry(:);
myTable.MaxJobs = repmat(maxjobs,height(myTable),1);
myTable.Node = fullname(:);
myTable.Object = object(:);
self.lqnGraph = self.lqnGraph.addnode(myTable);
self.nodeDep = nodeDep;
proc = self.hosts;
EndNodes = [];
Weight = [];
EdgeType = [];
PostType = [];
PreType = [];
for p=1:length(proc)
    tasks_p = proc{p}.tasks;
    for t=1:length(tasks_p)
        EndNodes(end+1,1) = findstring(self.lqnGraph.Nodes.Node,proc{p}.tasks(t).name);
        EndNodes(end,2) = findstring(self.lqnGraph.Nodes.Node,proc{p}.name);
        Weight(end+1,1) = 1.0;
        EdgeType(end+1,1) = 0; % within task
        PreType(end+1,1) = 0; % pre
        PostType(end+1,1) = 0; % post
        
        entries_tp = tasks_p(t).entries;
        for e=1:length(entries_tp)
            EndNodes(end+1,1) = findstring(self.lqnGraph.Nodes.Node,tasks_p(t).name);
            EndNodes(end,2) = findstring(self.lqnGraph.Nodes.Node,entries_tp(e).name);
            Weight(end+1,1) = 1.0;
            EdgeType(end+1,1) = 0; % within task
            PreType(end+1,1) = 0; % pre
            PostType(end+1,1) = 0; % post
        end
        
        sw_act_tp = tasks_p(t).activities; % sw activities
        for a=1:length(sw_act_tp)
            if ~isempty(sw_act_tp(a).boundToEntry)
                EndNodes(end+1,1) = findstring(self.lqnGraph.Nodes.Node,sw_act_tp(a).boundToEntry);
                EndNodes(end,2) = findstring(self.lqnGraph.Nodes.Node,sw_act_tp(a).name);
                Weight(end+1,1) = 1.0;
                EdgeType(end+1,1) = 0; % within task
                PreType(end+1,1) = 0; % pre
                PostType(end+1,1) = 0; % post
            end
            
            syncCallDests = sw_act_tp(a).syncCallDests;
            for sd=1:length(syncCallDests)
                EndNodes(end+1,1) = findstring(self.lqnGraph.Nodes.Node,sw_act_tp(a).name);
                EndNodes(end,2) = findstring(self.lqnGraph.Nodes.Node,syncCallDests{sd});
                Weight(end+1,1) = sw_act_tp(a).syncCallMeans(sd);
                EdgeType(end+1,1) = 1; % sync
                PreType(end+1,1) = 0; % pre
                PostType(end+1,1) = 0; % post
            end
            
            asyncCallDests = sw_act_tp(a).asyncCallDests;
            for asd=1:length(asyncCallDests)
                EndNodes(end+1,1) = findstring(self.lqnGraph.Nodes.Node,sw_act_tp(a).name);
                EndNodes(end,2) = findstring(self.lqnGraph.Nodes.Node,asyncCallDests{asd});
                Weight(end+1,1) = sw_act_tp(a).asyncCallMeans(asd);
                EdgeType(end+1,1) = 2; % async
                PreType(end+1,1) = 0; % pre
                PostType(end+1,1) = 0; % post
            end
        end
        
        act_prec_tp = tasks_p(t).precedences;
        for ap=1:length(act_prec_tp)
            if strcmpi(act_prec_tp(ap).postType,ActivityPrecedence.POST_SEQ)
                switch act_prec_tp(ap).preType
                    case ActivityPrecedence.PRE_SEQ
                        EndNodes(end+1,1) = findstring(self.lqnGraph.Nodes.Node, act_prec_tp(ap).preActs{1});
                        EndNodes(end,2) = findstring(self.lqnGraph.Nodes.Node, act_prec_tp(ap).postActs{1});
                        Weight(end+1,1) = 1.0;
                        EdgeType(end+1,1) = 0; % within task
                        PreType(end+1,1) = 0; % pre
                        PostType(end+1,1) = 0; % post
                    case ActivityPrecedence.PRE_AND
                        for pra=1:length(act_prec_tp(ap).preActs)
                            EndNodes(end+1,1) = findstring(self.lqnGraph.Nodes.Node, act_prec_tp(ap).preActs{pra});
                            EndNodes(end,2) = findstring(self.lqnGraph.Nodes.Node, act_prec_tp(ap).postActs{1});
                            if ~isempty(act_prec_tp(ap).preParams)
                                Weight(end+1,1) = act_prec_tp(ap).preParams(1)/length(act_prec_tp(ap).preActs);
                            else
                                Weight(end+1,1) = 1.0;
                            end
                            EdgeType(end+1,1) = 0; % within task
                            PreType(end+1,1) = 1; % pre-AND
                            PostType(end+1,1) = 0; % post
                        end
                    case ActivityPrecedence.PRE_OR
                        for pra=1:length(act_prec_tp(ap).preActs)
                            EndNodes(end+1,1) = findstring(self.lqnGraph.Nodes.Node, act_prec_tp(ap).preActs{pra});
                            EndNodes(end,2) = findstring(self.lqnGraph.Nodes.Node, act_prec_tp(ap).postActs{1});
                            Weight(end+1,1) = 1.0;
                            EdgeType(end+1,1) = 0; % within task
                            PreType(end+1,1) = 2; % pre-OR
                            PostType(end+1,1) = 0; % post
                        end
                    otherwise
                        line_error(mfilename,'Precedence is not supported yet.');
                end
            elseif strcmpi(act_prec_tp(ap).preType,ActivityPrecedence.PRE_SEQ)
                switch act_prec_tp(ap).postType
                    case ActivityPrecedence.POST_AND
                        for poa=1:length(act_prec_tp(ap).postActs)
                            EndNodes(end+1,1) = findstring(self.lqnGraph.Nodes.Node, act_prec_tp(ap).preActs{1});
                            EndNodes(end,2) = findstring(self.lqnGraph.Nodes.Node, act_prec_tp(ap).postActs{poa});
                            Weight(end+1,1) = 1.0;
                            EdgeType(end+1,1) = 0; % within task
                            PreType(end+1,1) = 0; % pre
                            PostType(end+1,1) = 1; % post-AND
                        end
                    case ActivityPrecedence.POST_OR
                        for poa=1:length(act_prec_tp(ap).postActs)
                            EndNodes(end+1,1) = findstring(self.lqnGraph.Nodes.Node, act_prec_tp(ap).preActs{1});
                            EndNodes(end,2) = findstring(self.lqnGraph.Nodes.Node, act_prec_tp(ap).postActs{poa});
                            Weight(end+1,1) = act_prec_tp(ap).postParams(poa);
                            EdgeType(end+1,1) = 0; % within task
                            PreType(end+1,1) = 0; % pre
                            PostType(end+1,1) = 2; % post-OR
                        end
                    case ActivityPrecedence.POST_LOOP
                        for poa=1:length(act_prec_tp(ap).postActs)
                            EndNodes(end+1,1) = findstring(self.lqnGraph.Nodes.Node, act_prec_tp(ap).preActs{1});
                            EndNodes(end,2) = findstring(self.lqnGraph.Nodes.Node, act_prec_tp(ap).postActs{poa});
                            if poa < length(act_prec_tp(ap).postActs)
                                Weight(end+1,1) = act_prec_tp(ap).postParams(poa);
                            else
                                Weight(end+1,1) = 1.0;
                            end
                            EdgeType(end+1,1) = 0; % within task
                            PreType(end+1,1) = 0; % pre
                            PostType(end+1,1) = 3; % post-LOOP
                        end
                    otherwise
                        
                        line_error(mfilename,'Precedence is not supported yet.');
                end
            else
                line_error(mfilename,'Precedence is not supported yet.');
            end
        end
    end
end

for e=1:size(EndNodes,1)
    source = EndNodes(e,1);
    target = EndNodes(e,2);
    self.lqnGraph = self.lqnGraph.addedge(name{source},name{target},Weight(e,1));
end

% do not merge with the previous loop as addedge does not preserve the order
for e=1:size(EndNodes,1)
    source = EndNodes(e,1);
    target = EndNodes(e,2);
    eid = self.lqnGraph.findedge(name{source},name{target});
    self.endNodes(eid,1) = EndNodes(e,1);
    self.endNodes(eid,2) = EndNodes(e,2);
    self.lqnGraph.Edges.Type(eid)=EdgeType(e,1);
    self.lqnGraph.Edges.Pre(eid)=PreType(e,1);
    self.lqnGraph.Edges.Post(eid)=PostType(e,1);
end

% fix all activities
for a=find(cellfun(@(c) strcmpi(c,'NaN'),self.lqnGraph.Nodes.Entry))'
    self.lqnGraph.Nodes.Entry{a} = self.findEntryOfActivity(self.lqnGraph.Nodes.Name{a});
end

% put hosts as target of hosted task
keep = ([findstring(self.lqnGraph.Nodes.Type,'H'); findstring(self.lqnGraph.Nodes.Type,'R'); findstring(self.lqnGraph.Nodes.Type,'T')]);
self.taskGraph = self.lqnGraph.subgraph(keep(keep>0)); % ignore missing types
% now H contains all tasks and edges to the hosts they run on
% we add external calls

actset = findstring(self.lqnGraph.Nodes.Type,'A')';
for s = actset(actset>0)
    source_activity = self.getNodeName(s);
    source_task = self.getNodeTask(source_activity);
    entryset = self.lqnGraph.successors(s)';
    for t = entryset(entryset>0)
        if strcmpi(self.getNodeType(t),'E')
            target_task = self.getNodeTask(t);
            edge_index = self.findEdgeIndex(s, t);
            edge_data = self.lqnGraph.Edges(edge_index, 2:end);
            self.taskGraph = self.taskGraph.addedge(source_task, target_task, edge_data);
        end
    end
end

self.param.Nodes.RespT=zeros(height(self.lqnGraph.Nodes),1);
self.param.Nodes.Util=zeros(height(self.lqnGraph.Nodes),1);
self.param.Nodes.Tput=zeros(height(self.lqnGraph.Nodes),1);

self.param.Edges.RespT=zeros(height(self.lqnGraph.Edges),1);
self.param.Edges.Tput=zeros(height(self.lqnGraph.Edges),1);

if isempty(self.indexes.reftasks)
    self.indexes.reftasks = findstring(self.lqnGraph.Nodes.Type, 'R');
end
if isempty(self.indexes.entries)
    self.indexes.entries = findstring(self.lqnGraph.Nodes.Type, 'E');
end
if isempty(self.indexes.hosts)
    self.indexes.hosts = findstring(self.lqnGraph.Nodes.Type, 'H');
end
if isempty(self.indexes.tasks)
    self.indexes.tasks = findstring(self.lqnGraph.Nodes.Type, 'T');
end
if isempty(self.indexes.activities)
    self.indexes.activities = findstring(self.lqnGraph.Nodes.Type, 'A');
end

end
