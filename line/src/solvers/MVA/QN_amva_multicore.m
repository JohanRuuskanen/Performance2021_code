function r = QN_amva_multicore(n,c)
% R = QN_AMVA_MULTICORE(N,C)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

M = length(n);
r = zeros(1,M);

for i = 1:M
    if isinf(c(i)) % delay server
        r(i) = 1/n(i);
    else % regular server with c(i) servers
        alpha=20;
        r(i) = 1./approximate(n(i),c(i),alpha);
        if isnan(r(i)) % if numerical problems
            r(i) = min(n(i),c(i));
        end
    end
end

end

function value = approximate(x,y,alpha) %soft-min
% VALUE = APPROXIMATE(X,Y,ALPHA) %SOFT-MIN
value = - ((-x)*exp(-alpha*x) -y*exp(-alpha*y)) / (exp(-alpha*x) + exp(-alpha*y));
end
