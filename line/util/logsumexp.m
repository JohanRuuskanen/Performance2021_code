function lse = logsumexp(x)
% L = LOGSUMEXP(X)
% Approximate the logarithm of a sum of exponentials
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

xstar = max(x);
lse = xstar + log(sum(exp(x-xstar)));
end