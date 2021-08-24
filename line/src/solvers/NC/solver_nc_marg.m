function [Pr,G,runtime] = solver_nc_marg(qn, options, lG)
% [PR,G,RUNTIME] = SOLVER_NC_MARG(QN, OPTIONS)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

if nargin == 2
    lG = NaN;
end

M = qn.nstations;    %number of stations
K = qn.nclasses;    %number of classes
state = qn.state;
S = qn.nservers;
NK = qn.njobs';  % initial population per class
C = qn.nchains;

PHr = qn.proc;

%% initialization

% determine service times
ST = zeros(M,K);
V = zeros(M,K);
for k = 1:K
    for i=1:M
        ST(i,k) = 1 ./ map_lambda(PHr{i,k});
    end
end
ST(isnan(ST))=0;
for c=1:qn.nchains
    V = V + qn.visits{c};
end

alpha = zeros(qn.nstations,qn.nclasses);
Vchain = zeros(qn.nstations,qn.nchains);
REFchain = zeros(qn.nstations,qn.nchains);
for c=1:qn.nchains
    inchain = find(qn.chains(c,:));
    for i=1:qn.nstations
        Vchain(i,c) = sum(qn.visits{c}(i,inchain)) / sum(qn.visits{c}(qn.refstat(inchain(1)),inchain));
        REFchain(i,c) = 1 / sum(qn.visits{c}(qn.refstat(inchain(1)),inchain));
        for k=inchain
            alpha(i,k) = alpha(i,k) + qn.visits{c}(i,k) / sum(qn.visits{c}(i,inchain)); % isn't alpha(i,j) always zero when entering here?
        end
    end
end

Vchain(~isfinite(Vchain))=0;
alpha(~isfinite(alpha))=0;

Lchain = zeros(M,C);
STchain = zeros(M,C);

Nchain = zeros(1,C);
refstatchain = zeros(C,1);
for c=1:qn.nchains
    inchain = find(qn.chains(c,:));
    isOpenChain = any(isinf(qn.njobs(inchain)));
    for i=1:qn.nstations
        % we assume that the visits in L(i,inchain) are equal to 1
        STchain(i,c) = ST(i,inchain) * alpha(i,inchain)';
        if isOpenChain && i == qn.refstat(inchain(1)) % if this is a source ST = 1 / arrival rates
            STchain(i,c) = 1 / sumfinite(qn.rates(i,inchain)); % ignore degenerate classes with zero arrival rates
        else
            STchain(i,c) = ST(i,inchain) * alpha(i,inchain)';
        end
    end
    Nchain(c) = sum(NK(inchain));
    refstatchain(c) = qn.refstat(inchain(1));
    if any((qn.refstat(inchain(1))-refstatchain(c))~=0)
        line_error(sprintf('Classes in chain %d have different reference station.',c));
    end
end
STchain(~isfinite(STchain))=0;
Tstart = tic;

[M,K]=size(STchain);

Lchain = zeros(M,K);
mu = ones(M,sum(Nchain));
for i=1:M
    Lchain(i,:) = STchain(i,:) .* Vchain(i,:);
    if isinf(S(i)) % infinite server
        mu(i,1:sum(Nchain)) = 1:sum(Nchain);
    else
        mu(i,1:sum(Nchain)) = min(1:sum(Nchain), S(i)*ones(1,sum(Nchain)));
    end
end
Lchain(~isfinite(Lchain))=0;

if isnan(lG)
    G = pfqn_gmvald(Lchain, Nchain, mu);
else
    G = exp(lG);
end

for ist=1:qn.nstations
    ind = qn.stationToNode(ist);
    isf = qn.stationToStateful(ist);
    [~,nirvec,sivec,kirvec] = State.toMarginal(qn, ind, state{isf});
    if min(nirvec) < 0 % user flags that state of i should be ignored
        Pr(i) = NaN;
    else
        set_ist = setdiff(1:qn.nstations,ist);
        nivec_chain = nirvec * qn.chains';
        G_minus_i = pfqn_gmvald(Lchain(set_ist,:), Nchain-nivec_chain, mu(set_ist,:), options);
        F_i = 1;
        switch qn.schedid(ist)
            case SchedStrategy.ID_FCFS
                for r=1:K
                    PHr = qn.proc{ist,r};
                    if ~isempty(PHr)
                        kir = kirvec(1,r,:); kir=kir(:)';
                        if length(kir)>1
                            line_error(mfilename,'Cannot return state probability because the product-form solution requires exponential service times at FCFS nodes.');
                        end
                        if ST(ist,r)~=max(ST(ist,:))
                            line_error(mfilename,'Cannot return state probability because the product-form solution requires identical service times at FCFS nodes.');
                        end
                    end
                end
                ci = find(sivec);
                if ~isempty(ci)
                    F_i = F_i * prod(exp(nirvec(1,:).*log(V(ist,r))))./prod(mu(ist,1:sum(kirvec(:))));
                else
                    F_i = 1;
                end
            case SchedStrategy.ID_SIRO
                for r=1:K
                    PHr = qn.proc{ist,r};
                    if ~isempty(PHr)
                        kir = kirvec(1,r,:); kir=kir(:)';
                        if length(kir)>1
                            line_error(mfilename,'Cannot return state probability because the product-form solution requires exponential service times at RAND nodes.');
                        end
                        if ST(ist,r)~=max(ST(ist,:))
                            line_error(mfilename,'Cannot return state probability because the product-form solution requires identical service times at RAND nodes.');
                        end
                    end
                end
                ci = find(sivec);
                if ~isempty(ci)
                    F_i = (nirvec(ci)/sum(nirvec)) * pfqn_gmvald(ST(ist,:).*V(ist,:), nirvec, mu(ist,:), options);
                else
                    F_i = 1;
                end
            case SchedStrategy.ID_PS
                for r=1:K
                    PHr = qn.proc{ist,r};
                    if ~isempty(PHr)
                        kir = kirvec(1,r,:); kir=kir(:)';
                        Ar = map_pie(PHr)*inv(-PHr{1});
                        F_i = F_i * prod(exp(kir.*log(V(ist,r)*Ar)))./prod(factorial(kir));
                    end
                end
                F_i = F_i * factorial(sum(kirvec(:)))./prod(mu(ist,1:sum(kirvec(:))));
            case SchedStrategy.ID_INF
                for r=1:K
                    PHr = qn.proc{ist,r};
                    if ~isempty(PHr)
                        kir = kirvec(1,r,:); kir=kir(:)';
                        Ar = map_pie(PHr)*inv(-PHr{1});
                        F_i = F_i * prod(exp(kir.*log(V(ist,r)*Ar)))./prod(factorial(kir));
                    end
                end
        end
        Pr(ist) =  F_i * G_minus_i / G;
    end
end

runtime = toc(Tstart);
Pr(isnan(Pr))=0;
lG = log(G);
if options.verbose > 0
    line_printf('\nNormalizing constant (NC) analysis completed. Runtime: %f seconds.\n',runtime);
end
return
end
