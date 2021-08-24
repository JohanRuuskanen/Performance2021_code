function [XN,QN,UN,CN,lGN,isNumStable,pi]=pfqn_mvald(L,N,Z,mu)
% [XN,QN,UN,CN,LGN]=PFQN_MVALD(L,N,Z,MU)

warn = true;
isNumStable = true;
[M,R]=size(L); % get number of queues (M) and classes (R)
if ~exist('mi','var')
    mi = ones(M,1);
end
Xs=zeros(R,prod(N+1)); % throughput for a model with station i less
pi=ones(M,sum(N)+1,prod(N+1)); % marginal queue-length probabilities pi(k)
WN=zeros(M,R);
n=pprod(N); % initialize the current population
lGN = 0;
while n~=-1
    WN=0*WN;
    for s=1:R
        if n(s)>0
            for i=1:M
                WN(i,s)=0;
                for k=1:sum(n)
                    WN(i,s)=WN(i,s)+(L(i,s)/mu(i,k))*k*pi(i,(k-1)+1,hashpop(oner(n,s),N));
                end
            end
            Xs(s,hashpop(n,N))=n(s)/(Z(s)+sum(WN(:,s)));
        end
    end
    
    % compute pi(k|n)
    for k=1:sum(n)
        for i=1:M
            pi(i,(k)+1,hashpop(n,N))=0;
        end
        for s=1:R
            if n(s)>0
                for i=1:M
                    pi(i,(k)+1,hashpop(n,N)) = pi(i,(k)+1,hashpop(n,N)) + (L(i,s)/mu(i,k))*Xs(s,hashpop(n,N))*pi(i,(k-1)+1,hashpop(oner(n,s),N));
                end
            end
        end
    end
    
    % compute pi(0|n)
    for i=1:M
        p0 = 1-sum(pi(i,(1:sum(n))+1,hashpop(n,N)));
        if p0<eps 
            if warn
                line_warning(mfilename,'MVA-LD is numerically unstable on this model, forcing all probabilities to be non-negative.'); 
%                N
                warn=false;
                isNumStable = false;
            end
            pi(i,(0)+1,hashpop(n,N)) = eps;
        else
            pi(i,(0)+1,hashpop(n,N)) = p0;
        end
    end
    
    last_nnz = max(find(n>0));
    if sum(n(1:last_nnz-1)) == sum(N(1:last_nnz-1)) & sum(n((last_nnz+1):R))==0
        logX = log(Xs(last_nnz,hashpop(n,N)));
        %hashpop(n,N)
        if ~isempty(logX)
            lGN(end+1) = lGN(end) - logX;
        end
    end
    
    n=pprod(n,N); % get the next population
end
X=Xs(:,hashpop(N,N));
XN=X';
pi=pi(:,:,hashpop(N,N));
QN = WN.*repmat(XN,M,1);
%UN = repmat(XN,M,1) .* L;
UN = 1-pi(:,1);
CN = N./XN - Z; % cycle time exclusive of think time
end

function i=hashpop(n,N)
% I=HASHPOP(N,N)

i=1; % index of the empty population
for r=1:length(N)
    i=i+prod(N(1:r-1)+1)*n(r);
end
end
