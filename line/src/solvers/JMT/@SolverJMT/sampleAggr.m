function stationStateAggr = sampleAggr(self, node, numEvents)
% STATIONSTATEAGGR = SAMPLEAGGR(NODE, NUMEVENTS)

if ~exist('node','var')
    line_error(mfilename,'sampleAggr requires to specify a node.');
end
if ~exist('numEvents','var')
    numEvents = -1;
else
    line_warning(mfilename,'JMT does not allow to fix the number of events for individual nodes. The number of returned events may be inaccurate.');
    numEvents = numEvents - 1; % we include the initialization as an event
end

Q = self.model.getAvgQLenHandles();
% create a temp model
modelCopy = self.model.copy;
modelCopy.resetNetwork;

% determine the nodes to logs
isNodeClassLogged = false(modelCopy.getNumberOfNodes, modelCopy.getNumberOfClasses);
ind = self.model.getNodeIndex(node.getName);
for r=1:modelCopy.getNumberOfClasses
    isNodeClassLogged(ind,r) = true;
end
% apply logging to the copied model
Plinked = self.model.getLinkedRoutingMatrix();
isNodeLogged = max(isNodeClassLogged,[],2);
logpath = tempdir;
modelCopy.linkAndLog(Plinked, isNodeLogged, logpath);
% simulate the model copy and retrieve log data
solverjmt = SolverJMT(modelCopy, self.getOptions);
if numEvents > 0
    solverjmt.maxEvents = numEvents*self.model.getNumberOfNodes*self.model.getNumberOfClasses;
else
    solverjmt.maxEvents = -1;
    numEvents = self.getOptions.samples;
end
solverjmt.getAvg(); % log data
logData = SolverJMT.parseLogs(modelCopy, isNodeLogged, Metric.QLen);

% from here convert from nodes in logData to stations
qn = modelCopy.getStruct;
ind = self.model.getNodeIndex(node.getName);
isf = qn.nodeToStateful(ind);
t = [];
nir = cell(1,qn.nclasses);
event = cell(1,qn.nclasses);
%ids = cell(1,qn.nclasses);

for r=1:qn.nclasses
    if isempty(logData{ind,r})
        nir{r} = [];
    else
        [~,uniqTS] = unique(logData{ind,r}.t);
        if isNodeClassLogged(isf,r)
            if ~isempty(logData{ind,r})
                t = logData{ind,r}.t(uniqTS);
                t = [t(2:end);t(end)];
                nir{r} = logData{ind,r}.QLen(uniqTS);
                event{r} = logData{ind,r}.event;
                %ids{r} = logData{ind,r}.
            end
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

if length(t) < 1+numEvents
    line_warning(mfilename,'LINE could not estimate correctly the JMT simulation length to return the desired number of events at the specified node. Try to re-run increasing the number of events.');
end

stationStateAggr = struct();
stationStateAggr.handle = node;
stationStateAggr.t = t;
stationStateAggr.t = stationStateAggr.t(1:min(length(t),1+numEvents),:);
stationStateAggr.t = [0; stationStateAggr.t(1:end-1)];
stationStateAggr.state = cell2mat(nir);
stationStateAggr.state = stationStateAggr.state(1:min(length(t),1+numEvents),:);
%stationStateAggr.job_id = 

event = cellmerge(event);
event = {event{cellisa(event,'Event')}}';
event_t = cellfun(@(c) c.t, event);
event_t = event_t(event_t <= max(stationStateAggr.t));
[~,I]=sort(event_t);
stationStateAggr.event = {event{I}};
stationStateAggr.event = stationStateAggr.event';
stationStateAggr.isaggregate = true;
end