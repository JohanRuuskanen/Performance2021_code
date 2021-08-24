function retval = containsstr(varargin)
% R = CONTAINSSTR(STR)
% Determine if pattern is in string
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

if isoctave 
    retval = strfind(varargin{1},varargin{2});
else
    if length(varargin) == 2
        retval = builtin('contains',varargin{1},varargin{2});
    else
        retval = builtin('contains',varargin{1},varargin{2},varargin{3},varargin{4});
    end
end
end
