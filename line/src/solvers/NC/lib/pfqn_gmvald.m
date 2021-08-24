function [G,lG]=pfqn_gmvald(L,N,mu,options)
% [G,LG]=PFQN_GMVALD(L,N,MU,OPTIONS)

% G=pfqn_gmvald(L,N,mu)
% mu: MxN matrix of load-dependent rates
[M,R]=size(L);
if isempty(L)
    G = 0; lG = -Inf; return
end
if nargin==2
    mu=ones(M,sum(N));
end
if nargin<4
    options = SolverNC.defaultOptions;
end
isLoadDep = false;
isInfServer = [];
for i=1:M
    if min(mu(i,1:sum(N))) == 1 & max(mu(i,1:sum(N))) == 1
        isInfServer(i) = false;
        continue; % this is a LI station
    elseif all(mu(i,1:sum(N)) == 1:sum(N))
        isInfServer(i) = true;
        continue; % this is a infinite server station
    else
        isInfServer(i) = false;
        isLoadDep = true;
    end
end

if ~isLoadDep
    % if load-independent model then use faster pfqn_gmva solver
    Lli = L(find(~isInfServer),:);
    if isempty(Lli)
        Lli = 0*N;
    end
    Zli = L(find(isInfServer),:);
    if isempty(Zli)
        Zli = 0*N;
    end
    options.method='exact';
    G = pfqn_nc(Lli,N,Zli, options);
end

G=0;
if M==0 G=0; lG=log(G); return; end
if sum(N==zeros(1,R))==R G=1; lG=log(G); return; end

if R==1
    G=pfqn_gmvaldsingle(L,N,mu);
    lG=log(G);
    return
end

G=G + pfqn_gmvald(L(1:(M-1),:),N,mu(1:(M-1),:));
for r=1:R
    if N(r)>0
        if R>1
            N_1 = oner(N,r);
        else
            N_1 = N-1;
        end
        G = G + (L(M,r)/mu(M,1))*pfqn_gmvald(L,N_1,mushift(mu,M));
    end
end
lG=log(G);
return
end

function mushifted=mushift(mu,i)
% MUSHIFTED=MUSHIFT(MU,I)

% shifts the service rate vector
[M,N]=size(mu);

for m=1:M
    if m==i
        mushifted(m,1:(N-1))=mu(m,2:N);
    else
        mushifted(m,1:(N-1))=mu(m,1:(N-1));
    end
end
end

function G=pfqn_gmvaldsingle(L,N,mu)
% G=PFQN_GMVALDSINGLE(L,N,MU)

[M,R]=size(L);
if R>1
    line_error(mfilename,'multiclass model detected. gmvaldsingle is for single class models.');
end
g=L(1,1)*0;
for n=1:N
    g(0 +1,n +1, 1 +1)=0;
end
for m=1:M
    for tm=1:(N+1)
        g(m +1,0 +1,tm +1)=1;
    end
    for n=1:N
        for tm=1:(N-n+1)
            g(m +1, n +1, tm +1)= g(m-1 +1, n +1, 1 +1)+L(m)*g(m +1, n-1 +1, tm+1 +1)/mu(m,tm);
        end
    end
end
G=g(M +1,N +1,1 +1);
end
