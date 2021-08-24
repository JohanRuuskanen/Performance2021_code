function pos = matchrow(matrix, row)
% pos = matchrow(M, r)
% Return position of row r in matrix M if unique, -1 otherwise
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

if all(matrix(end,:) == row)
    pos = size(matrix,1);
else
    pos = find(all(bsxfun(@eq,matrix,row),2),1);
    if isempty(pos)
        pos = -1;
    end
end
end