function line_printf(MSG,varargin)
% LINE_PRINTF(MSG, VARARGIN)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

MSG = strrep(MSG,'\n','');
if ~contains(MSG,'...')
    if contains(MSG,'Summary')
        %cprintf('_black', sprintf('%s\n',sprintf(MSG, varargin{:})));
        fprintf(1, sprintf('%s\n',sprintf(MSG, varargin{:})));
    elseif contains(MSG,'Iter')
        %cprintf('_black', sprintf('%s',sprintf(MSG, varargin{:})));
        fprintf(1, sprintf('%s',sprintf(MSG, varargin{:})));
    else
        fprintf(1, '%s\n', sprintf(MSG, varargin{:}));
    end
else
    MSG = sprintf('%s',sprintf(MSG, varargin{:}));
    fprintf(1, '%s', MSG);
end
end
