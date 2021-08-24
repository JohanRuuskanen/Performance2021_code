function line_error(caller,errmsg)
% LINE_ERROR(CALLER, ERRMSG)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
error(sprintf('[%s] %s',caller,errmsg));
end
