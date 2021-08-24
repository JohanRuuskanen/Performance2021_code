function Pr = getProbAggr(self, node, state_a)
% PR = GETPROBSTATEAGGR(NODE, STATE_A)

if ~exist('state_a','var')
    state_a = self.model.getState{self.model.getStationIndex(node)};
end
stationStateAggr = self.sampleAggr(node);
rows = findrows(stationStateAggr.state, state_a);
t = stationStateAggr.t;
dt = [diff(t);0];
Pr = sum(dt(rows))/sum(dt);
end