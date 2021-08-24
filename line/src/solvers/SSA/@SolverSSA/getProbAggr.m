function ProbAggr = getProbAggr(self, node, state)
% PROBSTATEAGGR = GETPROBSTATEAGGR(NODE, STATE)

% we do not use probSysState as that is for joint states
TranSysStateAggr = self.sampleSysAggr;
isf = self.model.getStatefulNodeIndex(node);
TSS = cell2mat({TranSysStateAggr.t,TranSysStateAggr.state{isf}});
TSS(:,1)=[TSS(1,1);diff(TSS(:,1))];
if ~exist('state','var')
    state = self.model.getState{isf};
end
rows = findrows(TSS(:,2:end), state);
if ~isempty(rows)
    ProbAggr = sum(TSS(rows,1))/sum(TSS(:,1));
else
    line_warning(mfilename,'The state was not seen during the simulation.');
    ProbAggr = 0;
end
end