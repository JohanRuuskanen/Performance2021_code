function [Gn,lGn]=pfqn_gm(L,N,Z,order)
% [GN,LGN]=PFQN_GM(L,N,Z,ORDER)

% PFQN_GM Exact and approximate solution of closed product-form queueing
% networks by Grundmann-Moeller cubature rules
%
% [Gn,lGn]=pfqn_gm(L,N,Z,S)
% Input:
% L : MxR demand matrix. L(i,r) is the demand of class-r at queue i
% N : 1xR population vector. N(r) is the number of jobs in class r
% Z : 1xR think time vector. Z(r) is the total think time of class r
% S : degree of the cubature rule. Exact if S=ceil((sum(N)-1)/2).
%
% Output:
% Gn : estimated normalizing constat
% lGn: logarithm of Gn. If Gn exceeds the floating-point range, only lGn
%      will be correctly estimated.
%
% Reference:
% G. Casale. Accelerating performance inference over closed systems by
% asymptotic methods. ACM SIGMETRICS 2017.
% Available at: http://dl.acm.org/citation.cfm?id=3084445

[M,R]=size(L);

if isempty(L) || isempty(N) || sum(N)==0
    Gn=1;
    return;
end
if nargin<4
    order=ceil((sum(N)-1)/2);
end

if ~exist('Z','var')
    Nt=sum(N);
    beta=N/Nt;
    f=@(x) ([x(:);1-sum(x)]'*L);
    h=@(x) beta.*log(f(x));
    [I,Q,ns]=simplexquad(@(x) prod(exp(Nt*h(x))),M-1,order,1e-8);
    Gn=Q*exp(gammaln(1+sum(N)+M-1)-sum(gammaln(1+N)));
    Gn=Gn(end);
else
    steps = 1e4;
    Nt = sum(N);
    beta = N/Nt;
    Gn=0;
    vmax=Nt*10;
    tol = 1e-10; % break infinite integration if G doesn't grow beyond tol
    dv=vmax/steps;
    for v=0:dv:vmax
        Lv=L*v+repmat(Z,M,1);
        h = @(x) beta.*log(([x(:);1-sum(x)]'*Lv));
        [~,Q] = simplexquad(@(x) exp(sum(Nt*h(x))),M-1,order,1e-8);
        dG = exp(-v)*v^(M-1)*Q(end)*dv;
        Gn = Gn + dG;
        if v>0 && dG/Gn<tol
            break
        end
    end
    Gn = Gn * exp(-sum(factln(N)));
end
lGn = log(Gn);
end

function [I,Q,nv]=simplexquad(f,n,order,tol)
% [I,Q,NV]=SIMPLEXQUAD(F,N,ORDER,TOL)

[Q,nv]=grnmol(@(x)f(x),eye(n,n+1),order,tol);
I=Q(end);
end

function [Q,nv] = grnmol( f, V, s , tol)
% [Q,NV] = GRNMOL( F, V, S , TOL)

%
%   Q = grnmol( f, V )
%     computes approximation to the integral of f over an s-simplex
%     with vertices as columns of V, an n x (n+1) matrix, using
%     order 1, 2, ..., s (degree 2s+1) Grundmann-Moler rules.
%   Output Q is a vector approximations of degree 1, 3, ... 2s+1.
%   Example: % final two results should be 2/(11x10x9x8x7x6)
%    n = 4;  grnmol(@(x)x(1)^2*x(n)^5,eye(n,n+1),4) 
%     Reference:
%       "Invariant Integration Formulas for the N-Simplex by
%        Combinatorial Methods", A. Grundmann and H. M. Moller,
%        SIAM J Numer. Anal. 15(1978), pp. 282-290
%
[n,~] = size(V); % n is the dimension
d = 0; % order of the polynomial
Q = (zeros(1,s+1));
T = zeros(1,s);
Qv = Q;
T0 = tic;
Vol = 1/factorial((n));%abs(det(V(:,1:n)-V(:,np)*ones(1,n)))/prod(mp([1:n]));
nv = 0;
x=[];
while 1
    m = n + 2*d + 1;
    al = (ones(n,1));
    alz = 2*d + 1;
    Qs = 0;
    while 1
        x(end+1,:)=[V*[alz; al]/m;1-sum(V*[alz; al]/m)];
        Qs = Qs + f(V*[alz; al]/m);
        nv = nv + 1;
        for j = 1 : n
            alz = alz - 2;
            if alz > 0
                al(j) = al(j) + 2;
                break
            end
            alz = alz + al(j) + 1; al(j) = 1;
        end
        if alz == 2*d+1
            break
        end
    end
    d = d + 1;
    Qv(d) = Vol*Qs;
    Q(d) = (0);
    p = 2/(prod(2*[n+1:m]));
    for i = 1 : d
        Q(d) = Q(d) + ((m+2-2*i))^(2*d-1)*p*Qv(d+1-i);
        p = -p*(m+1-i)/i;
    end
    if d > s || (d>1 && abs(Q(d)-Q(d-1))<tol*Q(d-1))
        Q((d+1):end)=[];
        break,
    end
    T(d)=toc(T0);
end
end