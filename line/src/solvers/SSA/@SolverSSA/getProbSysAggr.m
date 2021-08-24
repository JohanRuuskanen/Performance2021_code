function ProbSysAggr = getProbSysAggr(self)
% PROBSYSSTATEAGGR = GETPROBSYSSTATEAGGR()

TranSysStateAggr = self.sampleSysAggr;
TSS = cell2mat([TranSysStateAggr.t,TranSysStateAggr.state(:)']);
TSS(:,1)=[TSS(1,1);diff(TSS(:,1))];
state = self.model.getState;
qn = self.model.getStruct;
nir = zeros(qn.nstateful,qn.nclasses);
for isf=1:qn.nstateful
    ind = qn.statefulToNode(isf);
    if size(state{isf},1) > 1
        line_warning(mfilename,'Some states at node %d will be ignored. Please assign the node with a specific state.', ind);
    end
    [~,nir(isf,:)] = State.toMarginal(qn, ind, state{isf}(1,:));
end
nir = nir';
rows = findrows(TSS(:,2:end), nir(:)');
if ~isempty(rows)
    ProbSysAggr = sum(TSS(rows,1))/sum(TSS(:,1));
else
    line_warning(mfilename,'The state was not seen during the simulation.');
    ProbSysAggr = 0;
end
end