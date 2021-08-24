function  RD = pfqn_stdf(L,N,Z,S,fcfsNodes,rates,tset)
% sojourn time distribution function at multiserver FCFS nodes
% J. McKenna 1987 JACM
stabilityWarnIssued = false;
[M,R]=size(L);
T=length(tset);
mu = zeros(M,sum(N));
for k=1:M
    for n=1:sum(N)
        mu(k,n) = min(S(k),n);
    end
end

% when t==0 the hkc function is not well defined
tset(tset == 0) = Distrib.Zero;

for k=fcfsNodes(:)
    if range(rates(k,:))>Distrib.Zero
        line_error(mfilename,'The FCFS stations has distinct service rates, the model is invalid.');
    end
    hMAP = cell(1,1+sum(N));
    for n=0:sum(N)
        if n < S(k)
            hMAP{1+n} = map_exponential(1/rates(k,1));
        else
            hMAP{1+n} = map_sumind({map_exponential(1/rates(k,1)), map_erlang((n-S(k)+1)/(S(k)*rates(k,1)),n-S(k)+1)});            
        end
        hkc(1:T,1+n) = map_cdf(hMAP{1+n}, tset)';
    end
    
    for r=1:R
        if L(k,r) > Distrib.Zero
            Nr = oner(N,r);
            [~,~,~,~,lGr,isNumStable]  = pfqn_mvald(L,Nr,Z,mu);
            lGr = lGr(end);
            if ~isNumStable && ~stabilityWarnIssued
                stabilityWarnIssued = true;
                line_warning(mfilename,'The computation of the sojourn time distribution is numerically unstable');
            end
            RD{k,r} = zeros(length(tset),2);
            RD{k,r}(:,2) = tset(:);
            
            %% this is the original code in the paper
            %             Gkrt = zeros(T,1);
            %             nvec = pprod(Nr);
            %             while nvec >= 0
            %                 [~,~,~,~,lGk,isNumStable]  = pfqn_mvald(L(setdiff(1:M,k),:),Nr-nvec,Z,mu(setdiff(1:M,k),:));
            %                 lGk = lGk(end);
            %                 if ~isNumStable && ~stabilityWarnIssued
            %                     stabilityWarnIssued = true;
            %                     line_warning(mfilename,'The computation of the sojourn time distribution is numerically unstable');
            %                 end
            %                 for t=1:T
            %                     if sum(nvec)==0
            %                         Gkrt(t) = Gkrt(t) + hkc(t,1+0) * exp(lGk);
            %                     else
            %                         lFk = nvec*log(L(k,:))' +factln(sum(nvec)) - sum(factln(nvec)) - sum(log(mu(k,1:sum(nvec))));
            %                         Gkrt(t) = Gkrt(t) + hkc(t,1+sum(nvec)) * exp(lFk + lGk);
            %                     end
            %                 end % nvec
            %                 nvec = pprod(nvec, Nr);
            %             end
            %             lGkrt = log(Gkrt);
            %             RD{k,r}(1:T,1) = exp(lGkrt - lGr);
            
            %% this is faster as it uses the recursive form for LD models
            Hkrt = [];
            [~,~,~,~,lGk]  = pfqn_mvald(L(setdiff(1:M,k),:),Nr,Z,mu(setdiff(1:M,k),:));
            lGk = lGk(end);
            for t=1:T
                [t,T]
                gammat = mu;
                for m=1:sum(Nr)
                    gammat(k,m) = mu(k,m) * hkc(t,1+m-1) / hkc(t,1+m);
                end
                gammak = mushift(gammat,k); gammak=gammak(:,1:(sum(Nr)-1));
                Hkrt(t) = hkc(t,1+0) * exp(lGk); % nvec = 0
                for s=1:R % nvec >= 1s
                    if Nr(s)>0
                        gammak                        
                        lYks_t  = pfqn_rd(L,oner(Nr,s),Z,gammak);
                        %RD = lYks_t
                        %NRP 
                        %lYks_t = pfqn_nrp(L,oner(Nr,s),Z,gammak);
                        %[~,~,~,~,lYks_t]  = pfqn_mvald(L,oner(Nr,s),Z,gammak); lYks_t = lYks_t(end);
                        %if t==10
                        %    keyboard
                        %end
                        Hkrt(t) = Hkrt(t) + (L(k,s) * hkc(t,1+0) / gammat(k,1)) * exp(lYks_t); % nvec = 0
                    end
                end
            end
            Hkrt(isnan(Hkrt)) = Distrib.Zero;
            lHkrt = log(Hkrt);
            
            RD{k,r}(1:T,1) = min(1,exp(lHkrt - lGr));
        end
    end
end
end

function Fmx=Fm(m,x)
if m==1
    Fmx = 1-exp(-x);
else
    A = 0;
    for j=0:(m-1)
        A = A + x^j/factorial(j);
    end
    Fmx = 1-exp(-x)*A;
end
end

function mushifted = mushift(mu,i)
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