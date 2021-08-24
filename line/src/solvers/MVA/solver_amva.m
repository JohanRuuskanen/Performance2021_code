function [Q,U,R,T,C,X] = solver_amva(ST,V,N,nservers,SCV,options,sched,schedparam,refstat)
% [Q,U,R,T,C,X] = SOLVER_AMVA(ST,V,N,S,SCV,OPTIONS,SCHED,SCHEDPARAM,REFSTAT)
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

[M,K]=size(ST);
% queue-dependent functions to capture multi-server and delay stations
ga = @(n) QN_amva_multicore(n,nservers);
Q0 = options.init_sol;
tol = options.iter_tol;
if nargin < 5 || isempty(Q0)
    % balanced initialization
    Q = ones(M,K);
    Q = Q ./ repmat(sum(Q,1),size(Q,1),1) .* repmat(N,size(Q,1),1);
    Q(isinf(Q))=0; % open classes
    Q(refstat(isinf(N)))=Inf;
else
    Q=Q0;
end

if ~exist('sched','var')
    sched = cell(M,1);
    for i=1:M
        if isinf(nservers(i))
            sched(i) = SchedStrategy.INF;
        else
            sched(i) = SchedStrategy.PS; % default for non-inf servers is PS
        end
    end
end

extSET = find(sched==SchedStrategy.EXT);
infSET = find(sched==SchedStrategy.INF);
dpsSET = find(sched==SchedStrategy.DPS);
fcfsSET = find(sched==SchedStrategy.FCFS);
if K==1
    pfSET = find(sched==SchedStrategy.SIRO | sched==SchedStrategy.PS | sched==SchedStrategy.FCFS);
else
    pfSET = find(sched==SchedStrategy.SIRO | sched==SchedStrategy.PS);
end


if ~exist('tol','var')
    tol = 1e-6;
end
Nt = sum(N(isfinite(N)));
delta  = (Nt - 1) / Nt;
deltar = (N - 1) ./ N;
deltar(isinf(N)) = 1;

Q_1 = ones(M,K,1)*Inf;
U = zeros(M,K);
T = zeros(M,K);
C = zeros(1,K);
W = zeros(M,K);
Uhi = zeros(M,K);
X = 1./sum(ST,1);
ocl = [];
if any(isinf(N))
    ocl = find(isinf(N));
    for r=ocl % open classes
        X(r) = 1 ./ ST(refstat(r),r);
        deltar(r) = 1;
    end
end
rset = setdiff(1:K,find(N==0));

%% inner iteration
iter = 0;
while max(max(abs(Q-Q_1))) > tol
    iter = iter + 1;
    if iter > options.iter_max
        break;
    end
    Q_1 = Q;
    for k=1:M
        for r=rset
            Ak{r}(k) = 1 + delta * sum(Q(k,:));
            Akr(k,r) = 1 + deltar(r) * Q(k,r);
        end
    end
    
    %    b = be(Akr);
    for r=rset
        g = ga(Ak{r});
        sd = setdiff(rset,r); % change here to add class priorities
        
        for k=infSET(:)'
            W(k,r) = ST(k,r);
        end
        
        for k=pfSET(:)'
            % C(k,r) = L(k,r) * g(k) * b(k,r) * (1 + delta * sum(Q(k,:))); % old - class-dependence
            if nservers(k)>1
                if ismember(r,ocl)
                    W(k,r) = ST(k,r) * g(k) * (1 + (Q(k,r)+sum(Q(k,sd))));
                else
                    W(k,r) = ST(k,r) * g(k) * (1 + delta * (Q(k,r)+sum(Q(k,sd))));
                end
            else
                W(k,r) = ST(k,r) * (1 + deltar(r) * Q(k,r) + sum(Q(k,sd)));
            end
        end
        
        for k=dpsSET(:)'
            if nservers(k)>1
                line_error(mfilename,'Multi-server DPS not supported yet in AMVA solver.')
            else
                tss=Inf; % time-scale separation threshold, was 5 now unused
                Uhi(k,r) = sum(U(k,schedparam(k,:)>tss*schedparam(k,r)));
                W(k,r) = ST(k,r) / (1-Uhi(k,r)) * (1 + deltar(r) * Q(k,r));
                for s=setdiff(sd,r)
                    if schedparam(k,s)==schedparam(k,r)
                        W(k,r) = W(k,r) + ST(k,r) / (1-Uhi(k,r)) * Q(k,s);
                    elseif schedparam(k,s)/schedparam(k,r)<tss
                        W(k,r) = W(k,r) + ST(k,r) / (1-Uhi(k,r)) * Q(k,s)*schedparam(k,s)/schedparam(k,r);
                    elseif schedparam(k,s)/schedparam(k,r)>tss
                        % if there is time-scale separation, do nothing
                        % all is accounted for by 1/(1-Uhi)
                    end
                end
            end
        end
        
        for k=fcfsSET(:)'
            if ST(k,r) == 0
                W(k,r) = 0;
            else
                if nservers(k)>1
                    B = ((deltar .* X .* V(k,:) .* ST(k,:)) / nservers(k)); % note: this is in 0-1 as a utilization
                else
                    B = ones(1,K);
                end
                W(k,r) = ST(k,r)*(nservers(k)-1)/nservers(k); % multi-server correction with serial think time
                W(k,r) = W(k,r) + (ST(k,r)/nservers(k)) * (1 + SCV(k,r))/2; % high SCV
                W(k,r) = W(k,r) + (1/nservers(k)) * (ST(k,r) * deltar(r) * Q(k,r)*B(r) + ST(k,sd).*B(sd)*Q(k,sd)'); % FCFS approximation + reducing backlog proportionally to server utilizations; somewhat similar to Rolia-Sevcik -  method of layers - Sec IV.
            end
        end
        
        if sum(W(:,r)) == 0
            X(r) = 0;
        else
            if isinf(N(r))
                C(r) = V(:,r)'*W(:,r);
                % X(r) remains constant
            elseif N(r)==0
                X(r) = 0;
                C(r) = 0;
            else
                C(r) = V(:,r)'*W(:,r);
                X(r) = N(r) / C(r);
            end
        end
        
        for k=1:M
            Q(k,r) = X(r) * V(k,r) * W(k,r);
            T(k,r) = X(r) * V(k,r);
        end
    end
    
    for k=1:M
        for r=rset
            if isinf(nservers(k)) % infinite server
                U(k,r) = V(k,r)*ST(k,r)*X(r);
            else
                U(k,r) = V(k,r)*ST(k,r)*X(r)/nservers(k);
            end
        end
    end
    
end

for k=1:M
    for r=1:K
        if V(k,r)*ST(k,r)>0
            switch sched(k)
                case {SchedStrategy.FCFS,SchedStrategy.PS,SchedStrategy.DPS}
                    if sum(U(k,:))>1
                        U(k,r) = min(1,sum(U(k,:))) * V(k,r)*ST(k,r)*X(r) / ((V(k,:).*ST(k,:))*X(:));
                    end
            end
        end
    end
end

R = Q./T;
X(~isfinite(X))=0;
U(~isfinite(U))=0;
Q(~isfinite(Q))=0;
R(~isfinite(R))=0;

X(N==0)=0;
U(:,N==0)=0;
Q(:,N==0)=0;
R(:,N==0)=0;
T(:,N==0)=0;
W(:,N==0)=0;

end
