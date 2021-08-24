function [outhash, outrate, outprob] =  afterEventHashed(qn, ind, inhash, event, class)
% [OUTHASH, OUTRATE, OUTPROB] =  AFTEREVENTHASHED(QN, IND, INHASH, EVENT, CLASS)

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
isSimulation = false;
[outspace, outrate, outprob] =  State.afterEvent(qn, ind, inspace, event, class, isSimulation);
if isempty(outspace)
    outhash = -1;
    outrate = 0;
    return
else
    outhash = State.getHash(qn, ind, outspace);
end
end
