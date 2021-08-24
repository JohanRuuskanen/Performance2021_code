function tranSysState = sampleSys(self, numSamples)
% TRANSYSSTATE = SAMPLESYS(NUMSAMPLES)
options = self.getOptions;
options.force = true;
if isempty(self.result) || ~isfield(self.result,'infGen')
    self.runAnalysis();
end
[infGen, eventFilt] = self.getGenerator();
stateSpace = self.getStateSpace();
stateSpaceAggr = self.getStateSpaceAggr();

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
tranSysState = struct();
tranSysState.handle = self.model.getStatefulNodes';
tranSysState.t = cumsum([0,sjt(1:end-1)']');
for isf=1:length(initState)
    tranSysState.state{isf} = stateSpace(sts,(nst(isf):nst(isf+1)-1));
end

tranSysState.event = {};
for e = 1:length(event)    
    for a=1:length(qn.sync{event(e)}.active)
        tranSysState.event{end+1} = qn.sync{event(e)}.active{a};
        tranSysState.event{end}.t = tranSysState.t(e);
    end
    for p=1:length(qn.sync{event(e)}.passive)
        tranSysState.event{end+1} = qn.sync{event(e)}.passive{p};
        tranSysState.event{end}.t = tranSysState.t(e);
    end
end
tranSysState.isaggregate = false;

end