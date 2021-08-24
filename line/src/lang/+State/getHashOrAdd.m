function [hashid, qn] = getHashOrAdd(qn, ind, inspace)
% [HASHID, QN] = GETHASHORADD(QN, IND, INSPACE)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

if isempty(inspace)
    hashid = -1;
    return
end

% ind: node index
%ist = qn.nodeToStation(ind);
isf = qn.nodeToStateful(ind);

if isempty(qn.space{isf})
    line_error(mfilename,'Station state space is not initialized. Use setStateSpace method.\n');
end

% resize
if size(inspace,2) < size(qn.space{isf},2)
    inspace = [zeros(size(inspace,1),size(qn.space{isf},2)-size(inspace,2)), inspace];
elseif size(inspace,2) > size(qn.space{isf},2)
    qn.space{isf} = [zeros(size(qn.space{isf},1),size(inspace,2)-size(qn.space{isf},2)),qn.space{isf}];
end

hashid = matchrows(qn.space{isf},inspace);
%hashid=zeros(size(inspace,1),1);
for j=1:size(inspace,1)
    %    hashid(j,1) = matchrow(qn.space{isf},inspace(j,:));
    if hashid(j,1) <0
        qn.space{isf}(end+1,:) = inspace(j,:);
        hashid(j,1) = size(qn.space{isf},1);
    end
end
end
