function line_warning(caller,MSG, varargin)
% LINE_WARNING(CALLER, ERRMSG)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

MSG=sprintf(MSG, varargin{:});
w = warning('QUERY','ALL');
switch w(1).state
    case 'on'
        warning('[%s] %s',caller,MSG);
        %line_printf(sprintf('Warning [%s]: %s\n',caller,errmsg));
    case 'off'
end
end
