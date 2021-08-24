function sampleAggr = sampleAggr(self, node, numSamples)
% S = SAMPLEAGGR(NODE, NUMSAMPLES)
options = self.getOptions;
options.force = true;
if isempty(self.result) || ~isfield(self.result,'infGen')
    self.runAnalysis();
end
[infGen, eventFilt] = self.getGenerator();
stateSpace = self.getStateSpace();

initState = self.model.getState;
nst = cumsum([1,cellfun(@length,initState)']);
s0 = cell2mat(initState(:)');

% set initial state
pi0 = zeros(1,size(stateSpace,1));
pi0(matchrow(stateSpace,s0))=1;

% filter all CTMC events as a marked Markovian arrival process
D1 = cellsum(eventFilt);
D0 = infGen-D1;
MMAP = mmap_normalize([{D0},{D1},eventFilt(:)']);

% now sampel the MMAP
[sjt,event,~,~,sts] = mmap_sample(MMAP,numSamples, pi0);

qn = self.model.getStruct;
sampleAggr = struct();
sampleAggr.handle = node;
sampleAggr.t = cumsum([0,sjt(1:end-1)']');
ind = self.model.getNodeIndex(node);
isf = qn.nodeToStateful(ind);
sampleAggr.state = stateSpace(sts,(nst(isf):nst(isf+1)-1));
[~,sampleAggr.state] = State.toMarginal(qn,qn.statefulToNode(isf),sampleAggr.state);

sampleAggr.event = {};
%nodeEvent = false(length(event),1);
%nodeTS = zeros(length(event),1);
for e = 1:length(event)
    for a=1:length(qn.sync{event(e)}.active)
        sampleAggr.event{end+1} = qn.sync{event(e)}.active{a};
        sampleAggr.event{end}.t = sampleAggr.t(e);
%        if  qn.sync{event(e)}.active{a}.node == ind 
%            nodeEvent(e) = true;
%            nodeTS(e) = tranSysState.t(e);
%        end
    end
    for p=1:length(qn.sync{event(e)}.passive)
        sampleAggr.event{end+1} = qn.sync{event(e)}.passive{p};
        sampleAggr.event{end}.t = sampleAggr.t(e);
%        if  qn.sync{event(e)}.passive{p}.node == ind 
%            nodeEvent(e) = true;
%            nodeTS(e) = tranSysState.t(e);
%        end
    end
end

%tranSysState.state = tranSysState.state([1;find(nodeTS>0)],:);
%tranSysState.t = unique(nodeTS);
%tranSysState.event = tranSysState.event(nodeEvent)';
sampleAggr.isaggregate = true;

end