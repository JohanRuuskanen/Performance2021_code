function self = disableMetric(self, Y)
% SELF = DISABLEMETRIC(Y)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

if iscell(Y)
    Y={Y{:}};
    for i=1:length(Y)
        Y{i}.disable();
    end
else
    Y.disable();
end
end
