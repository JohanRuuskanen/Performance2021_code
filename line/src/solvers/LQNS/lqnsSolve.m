function lqnsSolve(filename, maxIter, wantExact)
% LQNSSOLVE(FILENAME, MAXITER, WANTEXACT)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
if ~exist('wantExact','var')
    wantExact = false;
end
if ~exist('maxIter','var') || isempty(maxIter)
    if wantExact
        system(['lqns -Playering=srvn -Pmva=exact -x ',filename]);
    else
        system(['lqns -Playering=srvn -x ',filename]);
    end
else
    if wantExact
        system(['lqns -Playering=srvn -i',num2str(maxIter),' -Pmva=exact -x ',filename]);
    else
        system(['lqns -Playering=srvn -i',num2str(maxIter),' -x ',filename]);
    end
end
end
