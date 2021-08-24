function [XN,QN,UN,RN,numIters,AN]=pfqn_aql(L,N,Z,TOL,QN0)
[K,C]=size(L);
if nargin<3
    Z=zeros(1,K);
end
MAXITER=1e3;
if nargin<4
    TOL = 1e-7;
end

Q=cell(1,C+1);
R=cell(1,C+1);
X=cell(1,C+1);
gamma = zeros(K,C);

if nargin<5
    for t=0:C
        n=oner(N,t);
        for k=1:K
            Q{t+1}(k,1)=sum(n)/K;
        end
    end
else
    for t=0:C
        n=oner(N,t);
        for k=1:K
            Q{t+1}(k,1)=QN0(k);
        end
    end
    
end
it = 0;
while 1
    Q_olditer = Q;
    it = it + 1;
    for t=0:C
        n=oner(N,t);
        for k=1:K
            for s=1:C
                R{t+1}(k,s) = L(k,s)*(1+(sum(n)-1)*(Q{t+1}(k)/sum(n)-gamma(k,s)));
            end
        end
        for s=1:C
            X{t+1}(s) = n(s)/(Z(s)+sum(R{t+1}(:,s)));
        end
        for k=1:K
            Q{t+1}(k) = X{t+1}(:)'*R{t+1}(k,:)';
        end
    end % for t
    for k=1:K
        for s=1:C
            gamma(k,s) = (Q{0+1}(k)/sum(N)) - (Q{s+1}(k)/(sum(N)-1));
        end
    end
    
    if max(abs((Q_olditer{1}(:)-Q{1}(:))./Q{1}(:))) < TOL || it == MAXITER
        numIters=it;
        break
    end
end
XN = X{1};
RN = R{1};
for k=1:K
    for s=1:C
        UN(k,s) = XN(s)*L(k,s);
        QN(k,s) = UN(k,s)*(1+Q{s+1}(k));
        AN(k,s) = Q{s+1}(k);
    end
end
end