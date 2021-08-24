function stationStateAggr = sampleAggr(self, node, numSamples)
% SAMPLE = SAMPLEAGGR(NODE, NUMSAMPLES)
if ~exist('node','var')
    line_error(mfilename,'sampleAggr requires to specify a station.');
end

%if exist('numsamples','var')
    %line_warning(mfilename,'SolveSSA does not support the numsamples parameter, use instead the samples option upon instantiating the solver.');
%end

options = self.getOptions;
switch options.method
    case {'default','serial'}
        options.samples = numSamples;
        options.force = true;
        [~, tranSystemState, event] = self.runAnalysis(options);
        qn = self.model.getStruct;
        isf = self.model.getStatefulNodeIndex(node);
        [~,nir]=State.toMarginal(qn,qn.statefulToNode(isf),tranSystemState{1+isf});
        stationStateAggr = struct();
        stationStateAggr.handle = node;
        stationStateAggr.t = tranSystemState{1};
        stationStateAggr.state = nir;
        qn = self.model.getStruct;
        stationStateAggr.event = {};
        for e = 1:length(event)
            for a=1:length(qn.sync{event(e)}.active)
                stationStateAggr.event{end+1} = qn.sync{event(e)}.active{a};
                stationStateAggr.event{end}.t = stationStateAggr.t(e);
            end
            for p=1:length(qn.sync{event(e)}.passive)
                stationStateAggr.event{end+1} = qn.sync{event(e)}.passive{p};
                stationStateAggr.event{end}.t = stationStateAggr.t(e);
            end
        end       
        stationStateAggr.isaggregate = true;
    otherwise
        line_error(mfilename,'sampleAggr is not available in SolverSSA with the chosen method.');
end
stationStateAggr.t = [0; stationStateAggr.t(2:end)];
end