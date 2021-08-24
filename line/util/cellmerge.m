function d = cellmerge(c)
% D=CELLMERGE(C)
% Vertically stacks the matrix elements of a cell array C.
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
    d = cell(0);
    for i=1:length(c)
        if iscell(c{i})
            d(end+1:end+length(c{i}),:)=c{i};
        end
    end
end