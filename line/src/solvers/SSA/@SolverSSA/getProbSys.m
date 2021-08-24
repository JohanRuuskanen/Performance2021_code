function ProbSys = getProbSys(self)
% PROBSYSSTATE = GETPROBSYSSTATE()

TranSysState = self.sampleSys;
TSS = cell2mat([TranSysState.t,TranSysState.state(:)']);
TSS(:,1)=[TSS(1,1);diff(TSS(:,1))];
state = cell2mat(self.model.getState');
rows = findrows(TSS(:,2:end), state);
if ~isempty(rows)
    ProbSys = sum(TSS(rows,1))/sum(TSS(:,1));
else
    line_warning(mfilename,'The state was not seen during the simulation.');
    ProbSys = 0;
end
end