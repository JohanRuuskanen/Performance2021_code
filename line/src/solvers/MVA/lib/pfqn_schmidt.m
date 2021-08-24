function [XN,QN,UN,CN] = pfqn_schmidt(D,N,S,sched)
% [XN,QN,UN,CN] = PFQN_SCHMIDT(D,N,S,SCHED)

% utilization in general ld case does not work
[M,R] = size(D);
closedClasses = 1:R;
XN = zeros(1,R);
UN = zeros(M,R);
CN = zeros(M,R);
QN = zeros(M,R);
C = length(closedClasses); % number of closed classes
Dc = D(:,closedClasses);
Nc = N(closedClasses);
prods = zeros(1,C); % needed for fast hashing
for r=1:C
    prods(r)=prod(Nc(1:r-1)+1);
end
% Start at nc=(0,...,0)
kvec = pprod(Nc);
% Initialize L and Pc
L = {}; % mean queue-length
Pc = {}; % state probabilities
for i=1:M
    switch sched(i)
        case SchedStrategy.ID_INF
            L{i} = zeros(R, prod(1+Nc)); % mean queue-length
        case SchedStrategy.ID_PS
            if S(i) == 1
                L{i} = zeros(R, prod(1+Nc)); % mean queue-length
            else
                Pc{i} = zeros(1 + sum(Nc), prod(1+Nc)); % Pr(j|N)
            end
        case SchedStrategy.ID_FCFS
            if all(D(i,:)==D(i,1))
                if S(i) == 1
                    L{i} = zeros(R, prod(1+Nc)); % mean queue-length
                else
                    Pc{i} = zeros(1 + sum(Nc), prod(1+Nc)); % Pr(j|N)
                end
            else
                Pc{i} = zeros(prod(1+Nc), prod(1+Nc)); % Pr(jvec|N)
            end
    end
end
x = zeros(C,prod(1+Nc));
w = zeros(M,C,prod(1+Nc));
for i=1:M
    Pc{i}(1 + 0, hashpop(kvec,Nc,C,prods)) = 1.0;
end
u = zeros(M,C);
% Population recursion
while kvec>=0
    hkvec = hashpop(kvec,Nc,C,prods);
    nc = sum(kvec);
    kprods = zeros(1,C); % needed for fast hashing
    for r=1:C
        kprods(r)=prod(kvec(1:r-1)+1);
    end
    for i=1:M
        for c=1:C
            hkvec_c = hashpop(oner(kvec,c),Nc,C,prods);
            % Compute mean residence times
            for n=1:nc
                switch sched(i)
                    case SchedStrategy.ID_INF
                        w(i,c,hkvec) = D(i,c);
                    case SchedStrategy.ID_PS
                        if S(i) == 1
                            w(i,c,hkvec) = Dc(i,c) * (1 + L{i}(1 + n, hkvec_c));
                        else
                            w(i,c,hkvec) = (Dc(i,c) / S(i)) * (1 + L{i}(c,hkvec_c));
                            for i=0:S(i)-2
                                w(i,c,hkvec) = w(i,c,hkvec) + (S(i)-1-i)*Pc{i}(1+i, hkvec_c);
                            end
                        end
                    case SchedStrategy.ID_FCFS
                        if all(D(i,:)==D(i,1)) % product-form case
                            if S(i) == 1
                                w(i,c,hkvec) = Dc(i,c) * (1 + L{i}(1 + n, hkvec_c));
                            else
                                w(i,c,hkvec) = (Dc(i,c) / S(i)) * (1 + L{i}(c,hkvec_c));
                                for i=0:S(i)-2
                                    w(i,c,hkvec) = w(i,c,hkvec) + (S(i)-1-i)*Pc{i}(1+i, hkvec_c);
                                end
                            end
                        else
                            nvec = pprod(kvec);
                            while nvec >= 0
                                if sum(nvec) > 0
                                    hnvec_c = hashpop(oner(nvec,c),kvec,C,kprods);
                                    Bcn = D(i,c) + max(0,sum(nvec)-S(i))/(S(i)*(sum(nvec)-1)) * (nvec*D(i,:)' - D(i,c));
                                    w(i,c,hnvec) = w(i,c,hnvec) + Bcn * Pc{i}(hnvec_c, hkvec_c);
                                end
                                nvec = pprod(nvec, kvec);
                            end
                        end
                end
            end
        end
    end
    % Compute tput
    for c=1:C
        x(c,hkvec) = kvec(c) / sum(w(1:M,c,hkvec));
    end
    for i=1:M
        for c=1:C
            L{i}(c) = x(c,hkvec) * w(i,c,hkvec);
        end
        switch sched(i)
            case SchedStrategy.ID_PS
                if S(i) > 1
                    for n=1:min(S(i),sum(kvec))
                        for c=1:C
                            hkvec_c = hashpop(oner(kvec,c),Nc,C,prods);
                            Pc{i}(1 + n, hkvec) = Pc{i}(1 + n, hkvec) + Dc(i,c) * (1/n) * x(c,hkvec) * Pc{i}(1+(n-1), hkvec_c);
                        end
                    end
                    Pc{i}(1 + 0, hkvec) = max(eps,1-sum(Pc(i, 1 + (1:min(S(i),sum(kvec))), hkvec)));
                end
            case SchedStrategy.ID_FCFS
                if all(D(i,:)==D(i,1))
                    if S(i) > 1
                        for n=1:min(S(i),sum(kvec))
                            for c=1:C
                                hkvec_c = hashpop(oner(kvec,c),Nc,C,prods);
                                Pc{i}(1 + n, hkvec) = Pc{i}(1 + n, hkvec) + Dc(i,c) * (1/n) * x(c,hkvec) * Pc{i}(1+(n-1), hkvec_c);
                            end
                        end
                        if sum(kvec)>0
                            Pc{i}(1 + 0, hkvec) = max(eps,1-sum(Pc(i, 1 + (1:min(S(i),sum(kvec))), hkvec)));
                        end
                    end
                else
                    nvec = pprod(kvec);
                    while nvec >= 0
                        hnvec = hashpop(nvec,kvec,C,kprods);
                        if sum(nvec)>0
                            for c=1:C
                                hnvec_c = hashpop(oner(nvec,c),kvec,C,kprods);
                                hkvec_c = hashpop(oner(kvec,c),Nc,C,prods);
                                Bcn = D(i,c) + max(0,sum(nvec)-S(i))/(S(i)*(sum(nvec)-1)) * (nvec*D(i,:)' - D(i,c));
                                Pc{i}(hnvec, hkvec) = Pc{i}(hnvec, hkvec) + (1/nvec(c))*x(c,hkvec)*Bcn*Pc{i}(hnvec_c, hkvec_c);
                            end
                        end
                        nvec = pprod(nvec, kvec);
                    end
                end
        end
    end
    kvec = pprod(kvec, Nc);
end
end
