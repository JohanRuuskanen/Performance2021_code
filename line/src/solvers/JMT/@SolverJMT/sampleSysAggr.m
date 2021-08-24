function sysStateAggr = sampleSysAggr(self, numEvents)
% SYSSTATEAGGR = SAMPLESYSAGGR(NUMEVENTS)

if ~exist('numEvents','var')
    numEvents = self.options.samples;
end

numEvents = numEvents - 1; % we include the initialization as an event
Q = self.model.getAvgQLenHandles();
statStateAggr = cell(self.model.getNumberOfStations,1);
% create a temp model
modelCopy = self.model.copy;
modelCopy.resetNetwork;
qn = self.model.getStruct;

% determine the nodes to logs
isNodeClassLogged = false(modelCopy.getNumberOfNodes, modelCopy.getNumberOfClasses);
for i= 1:modelCopy.getNumberOfStations
    ind = self.model.getNodeIndex(modelCopy.getStationNames{i});
    if qn.nodetype(ind) ~= NodeType.Source
        for r=1:modelCopy.getNumberOfClasses
            if ~Q{i,r}.disabled || nargin == 1
                isNodeClassLogged(ind,r) = true;
            else
                isNodeClassLogged(node,r) = true;
            end
        end
    end
end
% apply logging to the copied model
Plinked = self.model.getLinkedRoutingMatrix();
isNodeLogged = max(isNodeClassLogged,[],2);
logpath = tempdir;
modelCopy.linkAndLog(Plinked, isNodeLogged, logpath);

% simulate the model copy and retrieve log data
options = self.getOptions; options.samples = numEvents;
solverjmt = SolverJMT(modelCopy, options);
solverjmt.maxEvents = numEvents;
solverjmt.runAnalysis(); % log data
logData = SolverJMT.parseLogs(modelCopy, isNodeLogged, Metric.QLen);

% from here convert from nodes in logData to stations
event = {};
qn = modelCopy.getStruct;
for ist= 1:qn.nstations
    isf = qn.stationToStateful(ist);
    ind = qn.stationToNode(ist);
    t = [];
    nir = cell(1,qn.nclasses);
    event{isf,r} = {};
    if qn.nodetype(ind) == NodeType.Source
        nir{r} = [];
    else
        for r=1:qn.nclasses
            if ~isempty(logData{isf,r})
                [~,uniqTSi] = unique(logData{isf,r}.t);
                if isNodeClassLogged(isf,r)
                    if ~isempty(logData{isf,r})
                        t = logData{isf,r}.t(uniqTSi);
                        event{isf,r} = logData{isf,r}.event;
                        %t = [t(2:end);t(end)];
                        nir{r} = logData{isf,r}.QLen(uniqTSi);
                    end
                end
            else
                nir{r} = [];
                event{isf,r} = {};
            end
        end
    end
    if isfinite(self.options.timespan(2))
        stopAt = find(t>self.options.timespan(2),1,'first');
        if ~isempty(stopAt) && stopAt>1
            t = t(1:(stopAt-1));
            for r=1:length(nir)
                nir{r} = nir{r}(1:(stopAt-1));
            end
        end
    end
    statStateAggr{ist} = struct();
    statStateAggr{ist}.handle = self.model.stations{ist};
    statStateAggr{ist}.t = t;
    statStateAggr{ist}.state = cell2mat(nir);
    statStateAggr{ist}.event = {event{isf,:}};
    statStateAggr{ist}.isaggregate = true;
    sysStateAggr.arv_job_id = logData{2}.arvID;
    sysStateAggr.dep_job_id = logData{2}.depID;
end

tranSysStateAggr = cell(1,1+self.model.getNumberOfStations);

tranSysStateAggr{1} = []; % timestamps
for i=1:self.model.getNumberOfStations % stations
    if isempty(tranSysStateAggr{1})
        tranSysStateAggr{1} = statStateAggr{i}.t;
    else
        tumax = min(max(tranSysStateAggr{1}),max(statStateAggr{i}.t));
        tranSysStateAggr{1} = union(tranSysStateAggr{1}, statStateAggr{i}.t);
        tranSysStateAggr{1} = tranSysStateAggr{1}(tranSysStateAggr{1}<=tumax);
        tranSysStateAggr{1} = union(tranSysStateAggr{1}, statStateAggr{i}.t);
    end
end

for i=1:self.model.getNumberOfStations % stations
    ind = qn.stationToNode(i);
    tranSysStateAggr{1+i} = [];
    [~,uniqTSi] = unique(statStateAggr{i}.t);
    if qn.nodetype(ind) ~= NodeType.Source
        for j=1:self.model.getNumberOfClasses % classes
            % we floor the interpolation as we hold the last state
            if ~isempty(uniqTSi)
                Qijt = interp1(statStateAggr{i}.t(uniqTSi), statStateAggr{i}.state(uniqTSi,j), tranSysStateAggr{1},'previous');
                if isnan(Qijt(end))
                    Qijt(end)=Qijt(end-1);
                end
                tranSysStateAggr{1+i} = [tranSysStateAggr{1+i}, Qijt];
            else
                Qijt = NaN*ones(length(tranSysStateAggr{1}),1);
                tranSysStateAggr{1+i} = [tranSysStateAggr{1+i}, Qijt];
            end
        end
    else
        tranSysStateAggr{1+i} = [Inf];
    end
end

sysStateAggr = struct();
sysStateAggr.handle = self.model.stations';
sysStateAggr.t = tranSysStateAggr{1};
sysStateAggr.state = {tranSysStateAggr{2:end}};

% % % now put the events in the .event cell
% eventSysStateAggr = cell(self.model.getNumberOfStations, qn.nclasses); % timestamps
% for i=1:self.model.getNumberOfStations % stations
%     for r=1:qn.nclasses
%         if isempty(statStateAggr{i}.event)
%             eventSysStateAggr{i,r} = struct();
%         else
%             eventSysStateAggr{i,r} = struct();
%             for e=1:size(statStateAggr{i}.event{r},1)
%                 if ~isempty(statStateAggr{i}.event{r}{e})
%                     %statStateAggr{i}.event{r}{e}.t
%                     %etpos = find(sysStateAggr.t == statStateAggr{i}.event{r}{e}.t);
%                     %if ~isempty(etpos)
%                     evtype = statStateAggr{i}.event{r}{e}.event;
%                     evfield = EventType.toText(evtype);
%                     if ~isfield(eventSysStateAggr{i,r},evfield)
%                         eventSysStateAggr{i,r}.(EventType.toText(evtype)) = {};
%                     end
%                     eventSysStateAggr{i,r}.(EventType.toText(evtype)){end+1,1} = statStateAggr{i}.event{r}{e};
%                     %end
%                 end
%             end
%             %eventSysStateAggr{i,r} = statStateAggr{i}.event{r};
%         end
%     end
% end
% %sysStateAggr.t = [0; sysStateAggr.t(1:end-1)];
% sysStateAggr.event = eventSysStateAggr;

event = cellmerge(event);
event = {event{cellisa(event,'Event')}}';
event_t = cellfun(@(c) c.t, event);
[~,I]=sort(event_t);
sysStateAggr.event = {event{I}};
sysStateAggr.event = sysStateAggr.event';
sysStateAggr.isaggregate = true;
sysStateAggr.arv_job_id = logData{2}.arvID;
sysStateAggr.dep_job_id = logData{2}.depID;
end