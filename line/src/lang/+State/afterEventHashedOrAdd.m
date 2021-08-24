function [outhash, outrate, outprob, qn] =  afterEventHashedOrAdd(qn, ind, inhash, event, class)
% [OUTHASH, OUTRATE, OUTPROB, QN] =  AFTEREVENTHASHEDORADD(QN, IND, INHASH, EVENT, CLASS)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
if inhash == 0
    outhash = -1;
    outrate = 0;
    return
end
% ind: node index
%ist = qn.nodeToStation(ind);
isf = qn.nodeToStateful(ind);

inspace = qn.space{isf}(inhash,:);
isSimulation = true; % allow state vector to grow, e.g. for FCFS buffers
[outspace, outrate, outprob] =  State.afterEvent(qn, ind, inspace, event, class,isSimulation);
if isempty(outspace)
    outhash = -1;
    outrate = 0;
    return
else
    [outhash, qn] = State.getHashOrAdd(qn, ind, outspace);
end
end
