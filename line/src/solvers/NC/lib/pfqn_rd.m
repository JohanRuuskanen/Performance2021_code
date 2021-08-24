function lGN = pfqn_rd(L,N,Z,mu)
[M,R]=size(L);
gamma = ones(M,sum(N));
mu = mu(M,1:sum(N));
s = sum(N)*ones(1,M);
for i=1:M
    gamma(i,:) = mu(i,:)/mu(i,s(i));    
end

beta = ones(M,sum(N));
for i=1:M
    beta(i,1) = gamma(i,1) / (1-gamma(i,1)) ;
    for j=2:sum(N)
        beta(i,j) = (1-gamma(i,j-1)) * (gamma(i,j) / (1-gamma(i,j)));
    end
end
beta(isnan(beta))=Inf;
beta(isinf(beta))=Inf;

y = L;
for i=1:M
    y(i,:) = y(i,:) / (mu(i,end));
end

C=0;
sld = s(s>1);
vmax = min(sum(sld-1),sum(N));

Y = pfqn_aql(y,N,Z);
[~,~,~,~,lEN,isNumStable] = pfqn_mvald(y*Y',vmax,0,beta);
%[~,~,lEN] = gmvaldsingle(y*Y',vmax,beta);
for vtot=0:vmax
    EN = exp(lEN(vtot+1));
    C = C + ((sum(N)-max(0,max(vtot-1)))/sum(N)) * EN;
end
%[~,lGsigma] = pfqn_mci(y,N,Z,1e5);
lGN = pfqn_nc(y,N,Z);
%[~,lGsigma] = pfqn_ca(y,N,Z);
lGN = lGN + log(C);
end

function [G,g,lGN]=gmvaldsingle(L,N,mu)
[M,R]=size(L);
g=zeros(M+1,N+1,1+1);
g(1,1,1)=L(1,1)*0;
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
G=g(M+1,N+1,1+1);
lGN=log(g(M+1,:,1+1));
end