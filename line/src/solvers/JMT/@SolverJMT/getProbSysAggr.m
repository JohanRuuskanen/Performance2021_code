function probSysStateAggr = getProbSysAggr(self)
% PROBSYSSTATEAGGR = GETPROBSYSSTATEAGGR()

TranSysStateAggr = self.sampleSysAggr;
TSS = cell2mat([TranSysStateAggr.t,TranSysStateAggr.state(:)']);
TSS(:,1)=[diff(TSS(:,1));0];
state = self.model.getState;
qn = self.model.getStruct;
nir = zeros(qn.nstateful,qn.nclasses);
for isf=1:qn.nstateful
    ind = qn.statefulToNode(isf);
    [~,nir(isf,:)] = State.toMarginal(qn, ind, state{isf});
end
nir = nir';
rows = findrows(TSS(:,2:end), nir(:)');
if ~isempty(rows)
    probSysStateAggr = sum(TSS(rows,1))/sum(TSS(:,1));
else
    line_warning(mfilename,'The state was not seen during the simulation.');
    probSysStateAggr = 0;
end
end