function hashid = getHash(qn, ind, inspace)
% HASHID = GETHASH(QN, IND, INSPACE)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

if isempty(inspace)
    hashid = -1;
    return
end

% ind: node index
%ist = qn.nodeToStation(ind);
isf = qn.nodeToStateful(ind);

inspace = [zeros(size(inspace,1),size(qn.space{isf},2)-size(inspace,2)), inspace];
if isempty(qn.space{isf})
    line_error(mfilename,'Station state space is not initialized. Use setStateSpace method.\n');
end
hashid=zeros(size(inspace,1),1);
for j=1:size(inspace,1)
    hashid(j,1) = matchrow(qn.space{isf},inspace(j,:));
end
end
