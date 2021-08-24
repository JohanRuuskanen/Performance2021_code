function tranSysState = sampleSys(self, numSamples)
% TRANSYSSTATE = SAMPLESYS(NUMSAMPLES)
options = self.getOptions;
if exist('numSamples','var')
    options.samples = numSamples;
else
    numSamples = options.samples;
end
switch options.method
    case {'default','serial'}
        [~, tranSystemState, tranSync] = self.runAnalysis(options);
        tranSysState = struct();
        tranSysState.handle = self.model.getStatefulNodes';
        tranSysState.t = tranSystemState{1};
        tranSysState.state = {tranSystemState{2:end}};
        tranSysState.event = tranSync;
        event = tranSysState.event;
        
        for i=1:size(tranSysState.state,2)
            if size(tranSysState.state{i},1) > numSamples
                tranSysState.t = tranSystemState(1:numSamples);
                tranSysState.state = tranSysState.state{i}(1:numSamples,:);
            end
        end
        qn = self.model.getStruct;
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
        
    otherwise
        line_error(mfilename,'sampleSys is not available in SolverSSA with the chosen method.');
end


end