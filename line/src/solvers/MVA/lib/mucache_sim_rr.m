function [pij,simdata]=mucache_sim_rr(lambda,m,accost,samples)
%samples=1e5
%M={};
%for iter_tr=1:100
    iter_tr
    M{iter_tr}=[];;
    u=size(lambda,1); % number of users
    n=size(lambda,2); % number of items
    h=size(lambda,3)-1; % number of lists
    
    mt = sum(m); % total capacity
    list =  zeros(1,n-mt); % list(i) is the list of item i
    pos_in_list =  sparse(n,h); % plist(i,l) is the position of item i in cache l
    
    % initialize positions
    pos_in_list(1:(n-mt),1)=1:(n-mt); % first n-m outside the cache
    for j=1:h % the others in order
        pos_in_list((length(list)+1):(length(list)+m(j)),1+j) = 1:m(j);
        list = [list, j*ones(1,m(j))];
    end
    
    lam = zeros(1,u*n); % arrival rates for all items in current state
    for v=1:u
        for k=1:n
            lam((v-1)*n+k) = lambda(v,k,1+list(k));
        end
    end
    lamt = sum(lam); % holding rate
    
    t = 0; % simulation time
    pij = zeros(n,1+h); % item position probabilities
    lU = log(rand(samples,1)); % random number stream
    %simdata=zeros(numevents,1+n);
    acsum = cell(u,n,1+h); % cumulative sum of access cost to same or next level
    for v=1:u
        for i=1:n
            for listi=0:h
                acsum{v,i,1+listi} = cumsum(accost{v,i}(1+listi, (1+listi):end));
            end
        end
    end
    
    misst = [];
    simdata=[];
    for it=1:samples
        %    if mod(it,1e6)==0
        %       it
        %       pij ./ repmat(sum(pij,2),1,1+h)
        %    end
        %    dt = exprnd(1/lamt);
        
        %% determine time to next request, save probability of current state
        dt = - lU(it)/lamt; % exponential random number (1-U) is also uniform
        t = t + dt;
        for k=1:length(list)
            pij(k,1+list(k)) = pij(k,1+list(k)) + dt;
        end
        
        %% select firing transition
        y = 1 + sum(rand*lamt > cumsum(lam));
        i = 1 + mod(y-1, n); % requested item
        v = ceil(y / n); % stream issuing the request
        
        %% update state
        listi = list(i); % list of requested item i
        if listi == 0
            misst(end+1) = t;
        end
        posi = pos_in_list(i,1+listi); % current position of i
        if listi<h
            %% choose target list
            acvec = acsum{v,i,1+listi};
            rvec = rand*ones(1,size(acvec,2));
            listk = listi-1+find(rvec <= acvec,1,'first');
            %% swap position
            posk = randi(m(listk)); % pick position of victim in target list
            k = find(pos_in_list(:,1+listk) == posk,1,'first'); % resolve the id of victim item
            list(i) = listk; % update new list of i
            list(k) = listi; % update new list of k
            pos_in_list(i,1+listi) = 0; % delete position of i
            pos_in_list(i,1+list(i)) = posk; % save position of i
            pos_in_list(k,1+listk) = 0; % delete position of k
            pos_in_list(k,1+list(k)) = posi; % update position of k
            %% update arrival rates
            for s=1:u % for all streams
                lamt = lamt - lam((s-1)*n+i) - lam((s-1)*n+k);
                lam((s-1)*n+i) = lambda(s,i,1+list(i));
                lam((s-1)*n+k) = lambda(s,k,1+list(k));
                lamt = lamt + lam((s-1)*n+i) + lam((s-1)*n+k);
            end
        end
        %    simdata(it,:) = [t,list];
        M{iter_tr}(end+1,:) = [t,length(misst) / t];
    end
    pij = pij ./ repmat(t,1,1+h);
%end
%M
%keyboard
end