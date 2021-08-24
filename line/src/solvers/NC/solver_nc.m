function [Qchain,Uchain,Rchain,Tchain,Cchain,Xchain] = solver_nc(STchain, Vchain, Nchain, S, gamma, refstatchain, options)
% [QCHAIN,UCHAIN,RCHAIN,TCHAIN,CCHAIN,XCHAIN] = SOLVER_NC(STCHAIN, VCHAIN, NCHAIN, S, GAMMA, REFSTATCHAIN, OPTIONS)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

[M,K]=size(STchain);

L = zeros(M,K);
Z = [];
infServers = [];
for i=1:M
    if isinf(S(i)) % infinite server
        Z(end+1,:) = STchain(i,:) .* Vchain(i,:);
        infServers(end+1) = i;
    else
        L(i,:) = STchain(i,:) .* Vchain(i,:);
    end
end
Nt = sum(Nchain(isfinite(Nchain)));

G = pfqn_gmva(L,Nchain,Z);
for r=1:K
    Xchain(r) = pfqn_gmva(L,oner(Nchain,r),Z) / G;
    for i=1:M
        Qchain(i,r) = L(i,r) * pfqn_gmva(Ladd(L,i),oner(Nchain,r),Z) / G;
    end
end

Rchain = Qchain ./ repmat(Xchain,M,1) ./ Vchain;
Rchain(infServers,:) = Z ./ Vchain(infServers,:);
Tchain = repmat(Xchain,M,1) .* Vchain;
Uchain = Tchain .* L;
Cchain = Nchain ./ Xchain - Z;
Qchain(infServers,:) = Z.*Xchain;

Xchain(~isfinite(Xchain))=0;
Uchain(~isfinite(Uchain))=0;
Qchain(~isfinite(Qchain))=0;
Rchain(~isfinite(Rchain))=0;

Xchain(Nchain==0)=0;
Uchain(:,Nchain==0)=0;
Qchain(:,Nchain==0)=0;
Rchain(:,Nchain==0)=0;
Tchain(:,Nchain==0)=0;

end
