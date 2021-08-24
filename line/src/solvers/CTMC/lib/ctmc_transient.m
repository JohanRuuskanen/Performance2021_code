function [pi,t]=ctmc_transient(Q,pi0,t0,t1)
% [PI,T]=CTMC_TRANSIENT(Q,PI0,T0,T1)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
if nargin==2
    t1=pi0;
    t0=0;
    pi0=ones(1,length(Q));pi0=pi0/sum(pi0);
end
if nargin==3
    t1=t0;
    t0=0;
end
[t,pi]=ode23(@ctmc_transientode,[t0,t1],pi0);
%%[t,pi]=ode45(@ctmc_transientode,[t0,t1],pi0); % standard order 4-5
%[t,pi]=ode113(@ctmc_transientode,[t0,t1],pi0);

    function dpidt=ctmc_transientode(t,pi)
        % DPIDT=CTMC_TRANSIENTODE(T,PI)
        
        pi=pi(:)';
        dpidt=pi*Q;
        dpidt=dpidt(:);
    end

end
