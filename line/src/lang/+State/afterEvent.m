function [outspace, outrate, outprob] =  afterEvent(qn, ind, inspace, event, class, isSimulation)
% [OUTSPACE, OUTRATE, OUTPROB] =  AFTEREVENT(QN, IND, INSPACE, EVENT, CLASS, ISSIMULATION)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

% outprob: event probability if the action is passive

%if ~exist('isSimulation','var')
%    isSimulation = false;
%end
M = qn.nstations;
R = qn.nclasses;
S = qn.nservers;
phasessz = qn.phasessz;
phaseshift = qn.phaseshift;
pie = qn.pie;
outspace = [];
outrate = [];
outprob = 1;
% ind: node index
isf = qn.nodeToStateful(ind);

hasOnlyExp = false; % true if all service processes are exponential
if qn.isstation(ind)
    ist = qn.nodeToStation(ind);
    K = phasessz(ist,:);
    Ks = phaseshift(ist,:);
    if max(K)==1
        hasOnlyExp = true;
    end
    mu = qn.mu;
    phi = qn.phi;
    ph = qn.proc;
    capacity = qn.cap;
    classcap = qn.classcap;
    if K(class) == 0 % if this class is not accepted at the resource
        return
    end
    V = sum(qn.nvars(ind,:));
    space_var = inspace(:,(end-V+1):end); % local state variables
    space_srv = inspace(:,(end-sum(K)-V+1):(end-V)); % server state
    space_buf = inspace(:,1:(end-sum(K)-V)); % buffer state
elseif qn.isstateful(ind)
    V = sum(qn.nvars(ind,:));
    % in this case service is always immediate so sum(K)=1
    space_var = inspace(:,(end-V+1):end); % local state variables
    space_srv = inspace(:,(end-R-V+1):(end-V)); % server state
    space_buf = []; % buffer state
else % stateless node
    space_var = [];
    space_srv = [];
    space_buf = [];
end

%switch qn.nodetype(ind)
%    case {NodeType.Queue, NodeType.Delay, NodeType.Source}
if qn.isstation(ind)
    switch event
        case EventType.ARV %% passive
            % return if there is no space to accept the arrival
            [ni,nir] = State.toMarginalAggr(qn,ind,inspace,K,Ks,space_buf,space_srv,space_var);
            % otherwise check scheduling strategy
            pentry = pie{ist,class};
            outprob = [];
            outprob_k = [];
            for kentry = 1:K(class)
                space_var_k = space_var;
                space_srv_k = space_srv;
                space_buf_k = space_buf;
                switch qn.schedid(ist)
                    case SchedStrategy.ID_EXT % source, can receive any "virtual" arrival from the sink as long as it is from an open class
                        if isinf(qn.njobs(class))
                            outspace = inspace;
                            outrate = -1*zeros(size(outspace,1)); % passive action, rate is unspecified
                            outprob = ones(size(outspace,1));
                            break
                        end
                    case {SchedStrategy.ID_PS, SchedStrategy.ID_INF, SchedStrategy.ID_DPS, SchedStrategy.ID_GPS}
                        % job enters service immediately
                        space_srv_k(:,Ks(class)+kentry) = space_srv_k(:,Ks(class)+kentry) + 1;
                        outprob_k = pentry(kentry)*ones(size(space_srv_k,1));
                    case {SchedStrategy.ID_SIRO, SchedStrategy.ID_SEPT, SchedStrategy.ID_LEPT}
                        if ni<S(ist)
                            space_srv_k(:,Ks(class)+kentry) = space_srv_k(:,Ks(class)+kentry) + 1;
                            outprob_k = pentry(kentry)*ones(size(space_srv_k,1));
                        else
                            space_buf_k(:,class) = space_buf_k(:,class) + 1;
                            outprob_k = pentry(kentry)*ones(size(space_srv_k,1));
                        end
                    case {SchedStrategy.ID_FCFS, SchedStrategy.ID_HOL, SchedStrategy.ID_LCFS}
                        % find states with all servers busy - this
                        % needs not to be moved
                        all_busy_srv = find(sum(space_srv_k,2) >= S(ist));
                        
                        % find and modify states with an idle server
                        idle_srv = sum(space_srv_k,2) < S(ist);
                        space_srv_k(idle_srv, end-sum(K)+Ks(class)+kentry) = space_srv_k(idle_srv,end-sum(K)+Ks(class)+kentry) + 1; % job enters service
                        
                        if isSimulation
                            if ni < capacity(ist) && nir(class) < classcap(ist,class) % if there is room
                                if ~any(space_buf_k(:)==0) % but the buffer has no empty slots
                                    % append job slot
                                    space_buf_k = [zeros(size(space_buf_k,1),1),space_buf_k];
                                end
                            end
                        end
                        %get position of first empty slot
                        empty_slots = -1*ones(all_busy_srv,1);
                        if size(space_buf_k,2) == 0
                            empty_slots(all_busy_srv) = false;
                        elseif size(space_buf_k,2) == 1
                            empty_slots(all_busy_srv) = space_buf_k(all_busy_srv,:)==0;
                        else
                            empty_slots(all_busy_srv) = max(bsxfun(@times, space_buf_k(all_busy_srv,:)==0, [1:size(space_buf_k,2)]),[],2);
                        end
                        
                        % ignore states where the buffer has no empty slots
                        wbuf_empty = empty_slots>0;
                        if any(wbuf_empty)
                            space_srv_k = space_srv_k(wbuf_empty,:);
                            space_buf_k = space_buf_k(wbuf_empty,:);
                            empty_slots = empty_slots(wbuf_empty);
                            space_buf_k(sub2ind(size(space_buf_k),1:size(space_buf_k,1),empty_slots')) = class;
                            %outspace(all_busy_srv(wbuf_empty),:) = [space_buf, space_srv, space_var];
                        end
                        outprob_k = pentry(kentry)*ones(size(space_srv_k,1));
                end
                outspace_k = [space_buf_k, space_srv_k, space_var_k];
                % remove states where new arrival violates capacity or cutoff constraints
                en = classcap(ist,class)> nir(:,class) | capacity(ist)*ones(size(ni,1),1) > ni;
                outspace = [outspace; outspace_k(en,:)];
                
                outrate = [outrate; -1*ones(size(outspace_k,1))]; % passive action, rate is unspecified
                outprob = [outprob; outprob_k];
            end
            if isSimulation
                if size(outprob,1) > 1
                    tot_prob = sum(outprob);
                    cum_prob = cumsum(outprob) / tot_prob;
                    firing_ctr = 1 + max([0,find( rand > cum_prob' )]); % select action
                    outspace = outspace(firing_ctr,:);
                    outrate = -1;
                    outprob = 1;
                end
            end
        case EventType.DEP
            if any(any(space_srv(:,(Ks(class)+1):(Ks(class)+K(class))))) % something is busy
                if hasOnlyExp && (qn.schedid(ist) == SchedStrategy.ID_PS || qn.schedid(ist) == SchedStrategy.ID_DPS || qn.schedid(ist) == SchedStrategy.ID_GPS || qn.schedid(ist) == SchedStrategy.ID_INF)
                    nir = space_srv;
                    ni = sum(nir,2);
                    sir = nir;
                    kir = sir;
                else
                    [ni,nir,sir,kir] = State.toMarginal(qn,ind,inspace,K,Ks,space_buf,space_srv,space_var);
                end
                switch qn.routing(ind)
                    case RoutingStrategy.ID_RRB
                        idx = find(space_var(end) == qn.varsparam{ind}.outlinks);
                        if idx < length(qn.varsparam{ind}.outlinks)
                            space_var = qn.varsparam{ind}.outlinks(idx+1);
                        else
                            space_var = qn.varsparam{ind}.outlinks(1);
                        end
                end
                if sir(class)>0 % is a job of class is in service
                    outprob = [];
                    for k=1:K(class)
                        space_srv = inspace(:,(end-sum(K)-V+1):(end-V)); % server state
                        space_buf = inspace(:,1:(end-sum(K)-V)); % buffer state
                        rate = zeros(size(space_srv,1),1);
                        en =  space_srv(:,Ks(class)+k) > 0;
                        if any(en)
                            switch qn.schedid(ist)
                                case SchedStrategy.ID_EXT % source, can produce an arrival from phase-k as long as it is from an open class
                                    if isinf(qn.njobs(class))
                                        pentry = pie{ist,class};
                                        for kentry = 1:K(class)
                                            space_srv = inspace(:,(end-sum(K)-V+1):(end-V)); % server state
                                            space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                            space_srv(en,Ks(class)+kentry) = space_srv(en,Ks(class)+kentry) + 1; % new job
                                            outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                            outrate = [outrate; pentry(kentry)*mu{ist,class}(k)*phi{ist,class}(k)*ones(size(inspace(en,:),1),1)];
                                            outprob = [outprob; ones(size(space_buf(en,:),1),1)];
                                        end
                                    end
                                case SchedStrategy.ID_INF % move first job in service
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    rate(en) = mu{ist,class}(k)*(phi{ist,class}(k)).*kir(en,class,k); % assume active
                                    % if state is unchanged, still add with rate 0
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(end,:)];
                                    outrate = [outrate; rate(en,:)];
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                case SchedStrategy.ID_PS % move first job in service
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    rate(en) = mu{ist,class}(k)*(phi{ist,class}(k)).*(kir(en,class,k)./ni(en)).*min(ni(en),S(ist)); % assume active
                                    % if state is unchanged, still add with rate 0
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    outrate = [outrate; rate(en,:)];
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                case SchedStrategy.ID_DPS
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    if S(ist) > 1
                                        line_error(mfilename,'Multi-server DPS stations are not supported yet.');
                                    end
                                    % in GPS, the scheduling parameter are the weights
                                    w_i = qn.schedparam(ist,:);
                                    w_i = w_i / sum(w_i);
                                    rate(en) = mu{ist,class}(k)*(phi{ist,class}(k))*(kir(en,class,k)/nir(class))*w_i(class)*nir(class)./(sum(repmat(w_i,sum(en),1)*nir',2));
                                    % if state is unchanged, still add with rate 0
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    outrate = [outrate; rate(en,:)];
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                case SchedStrategy.ID_GPS
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    if S(ist) > 1
                                        line_error(mfilename,'Multi-server GPS stations are not supported yet.');
                                    end
                                    % in GPS, the scheduling parameter are the weights
                                    w_i = qn.schedparam(ist,:);
                                    w_i = w_i / sum(w_i);
                                    cir = min(nir,ones(size(nir)));
                                    rate = mu{ist,class}(k)*(phi{ist,class}(k))*(kir(en,class,k)/nir(class))*w_i(class)/(w_i*cir(:)); % assume active
                                    % if state is unchanged, still add with rate 0
                                    outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                    outrate = [outrate; rate(en,:)];
                                    outprob = [outprob; ones(size(rate(en,:),1),1)];
                                case SchedStrategy.ID_FCFS % move first job in service
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    rate(en) = mu{ist,class}(k)*(phi{ist,class}(k)).*kir(en,class,k); % assume active
                                    en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                    en_wobuf = ~en_wbuf;
                                    outspace = [outspace; space_buf(en_wobuf,:), space_srv(en_wobuf,:), space_var(en_wobuf,:)];
                                    outrate = [outrate; rate(en_wobuf,:)];
                                    outprob = [outprob; ones(size(rate(en_wobuf,:),1),1)];
                                    if any(en_wbuf)
                                        start_svc_class = space_buf(en_wbuf,end);
                                        if start_svc_class > 0
                                            space_buf(en_wbuf,:) = [zeros(sum(en_wbuf),1),space_buf(en_wbuf,1:end-1)];
                                            pentry_svc_class = pie{ist,start_svc_class};
                                            for kentry = 1:K(start_svc_class)
                                                space_srv(en_wbuf,Ks(start_svc_class)+kentry) = space_srv(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                                outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                                rate_k = rate;
                                                rate_k(en_wbuf,:) = rate(en_wbuf,:)*pentry_svc_class(kentry);
                                                outrate = [outrate; rate_k(en,:)];
                                                outprob = [outprob; ones(size(rate(en,:),1),1)];
                                                space_srv(en_wbuf,Ks(start_svc_class)+kentry) = space_srv(en_wbuf,Ks(start_svc_class)+kentry) - 1;
                                            end
                                        end
                                    end
                                    % if state is unchanged, still add with rate 0
                                case SchedStrategy.ID_HOL % FCFS priority
                                    rate(en) = mu{ist,class}(k)*(phi{ist,class}(k)).*kir(:,class,k); % assume active
                                    en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                    en_wobuf = ~en_wbuf;
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    priogroup = [0,qn.classprio];
                                    space_buf_groupg = arrayfun(@(x) priogroup(1+x), space_buf);
                                    start_classprio = max(space_buf_groupg(en_wbuf,:),[],2);
                                    isrowmax = space_buf_groupg == repmat(start_classprio, 1, size(space_buf_groupg,2));
                                    [~,rightmostMaxPosFlipped]=max(fliplr(isrowmax),[],2);
                                    rightmostMaxPos = size(isrowmax,2) - rightmostMaxPosFlipped + 1;
                                    start_svc_class = space_buf(en_wbuf, rightmostMaxPos);
                                    outspace = [outspace; space_buf(en_wobuf,:), space_srv(en_wobuf,:), space_var(en_wobuf,:)];
                                    outrate = [outrate; rate(en_wobuf,:)];
                                    outprob = [outprob; ones(size(rate(en_wobuf,:),1),1)];
                                    if start_svc_class > 0
                                        pentry_svc_class = pie{ist,start_svc_class};
                                        for kentry = 1:K(start_svc_class)
                                            space_srv_k = space_srv;
                                            space_buf_k = space_buf;
                                            space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) = space_srv_k(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                            for j=find(en_wbuf)'
                                                space_buf_k(j,:) = [0, space_buf_k(j,1:rightmostMaxPos(j)-1), space_buf_k(j,(rightmostMaxPos(j)+1):end)];
                                            end
                                            % if state is unchanged, still add with rate 0
                                            outspace = [outspace; space_buf_k(en_wbuf,:), space_srv_k(en_wbuf,:), space_var(en_wbuf,:)];
                                            rate_k = rate;
                                            rate_k(en_wbuf,:) = rate(en_wbuf,:) * pentry_svc_class(kentry);
                                            outrate = [outrate; rate_k(en_wbuf,:)];
                                            outprob = [outprob; ones(size(rate_k(en_wbuf,:),1),1)];
                                        end
                                    end
                                case SchedStrategy.ID_LCFS % move last job in service
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    rate(en) = mu{ist,class}(k)*(phi{ist,class}(k)).*kir(:,class,k); % assume active
                                    en_wbuf = en & ni>S(ist); %states with jobs in buffer
                                    [~, colfirstnnz] = max( space_buf(en_wbuf,:) ~=0, [], 2 ); % find first nnz column
                                    start_svc_class = space_buf(en_wbuf,colfirstnnz); % job entering service
                                    space_buf(en_wbuf,colfirstnnz)=0;
                                    if isempty(start_svc_class)
                                        outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                        outrate = [outrate; rate(en,:)];
                                        outprob = [outprob; ones(size(rate(en,:),1),1)];
                                        return
                                    end
                                    for kentry = 1:K(start_svc_class)
                                        pentry_svc_class = pie{ist,start_svc_class};
                                        space_srv(en_wbuf,Ks(start_svc_class)+kentry) = space_srv(en_wbuf,Ks(start_svc_class)+kentry) + 1;
                                        % if state is unchanged, still add with rate 0
                                        outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                        rate_k = rate;
                                        rate_k(en_wbuf,:) = rate(en_wbuf,:)*pentry_svc_class(kentry);
                                        outrate = [outrate; rate_k(en,:)];
                                        outprob = [outprob; ones(size(rate(en,:),1),1)];
                                        space_srv(en_wbuf,Ks(start_svc_class)+kentry) = space_srv(en_wbuf,Ks(start_svc_class)+kentry) - 1;
                                    end
                                case SchedStrategy.ID_SIRO
                                    rate = zeros(size(space_srv,1),1);
                                    rate(en) = mu{ist,class}(k)*(phi{ist,class}(k)).*kir(:,class,k); % this is for states not in en_buf
                                    space_srv = inspace(:,end-sum(K)+1:end); % server state
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    % first record departure in states where the buffer is empty
                                    en_wobuf = en & sum(space_buf(en,:),2) == 0;
                                    outspace = [outspace; space_buf(en_wobuf,:), space_srv(en_wobuf,:), space_var(en_wobuf,:)];
                                    outrate = [outrate; rate(en_wobuf,:)];
                                    outprob = [outprob; ones(size(rate(en_wobuf,:),1),1)];
                                    % let's go now to states where the buffer is non-empty
                                    for r=1:R % pick a job of a random class
                                        rate_r = rate;
                                        space_buf = inspace(:,1:(end-sum(K))); % buffer state
                                        en_wbuf = en & space_buf(en,r) > 0; % states where the buffer is non-empty
                                        space_buf(en_wbuf,r) = space_buf(en_wbuf,r) - 1; % remove from buffer
                                        space_srv_r = space_srv;
                                        pentry_svc_class = pie{ist,r};
                                        pick_prob = (nir(r)-sir(r)) / (ni-sum(sir));
                                        if pick_prob >= 0
                                            rate_r(en_wbuf,:) = rate_r(en_wbuf,:) * pick_prob;
                                        end
                                        for kentry=1:K(r)
                                            space_srv_r(en_wbuf,Ks(r)+kentry) = space_srv_r(en_wbuf,Ks(r)+kentry) + 1; % bring job in service
                                            outspace = [outspace; space_buf(en_wbuf,:), space_srv_r(en_wbuf,:), space_var(en_wbuf,:)];
                                            rate_k = rate_r;
                                            rate_k(en_wbuf,:) = rate_k(en_wbuf,:) * pentry_svc_class(kentry);
                                            outrate = [outrate; rate_k(en_wbuf,:)];
                                            outprob = [outprob; ones(size(rate(en_wbuf,:),1),1)];
                                            space_srv_r(en_wbuf,Ks(r)+kentry) = space_srv_r(en_wbuf,Ks(r)+kentry) - 1; % bring job in service
                                        end
                                    end
                                case {SchedStrategy.ID_SEPT,SchedStrategy.ID_LEPT} % move last job in service
                                    rate = zeros(size(space_srv,1),1);
                                    rate(en) = mu{ist,class}(k)*(phi{ist,class}(k)).*kir(:,class,k); % this is for states not in en_buf
                                    space_srv = inspace(:,end-sum(K)+1:end); % server state
                                    space_srv(en,Ks(class)+k) = space_srv(en,Ks(class)+k) - 1; % record departure
                                    space_buf = inspace(:,1:(end-sum(K))); % buffer state
                                    % in SEPT, the scheduling parameter is the priority order of the class means
                                    % en_wbuf: states where the buffer is non-empty
                                    % sept_class: class to pick in service
                                    [en_wbuf, first_class_inrow] = max(space_buf(:,qn.schedparam(ist,:))~=0, [], 2);
                                    sept_class = qn.schedparam(ist,first_class_inrow); % this is different for sept and lept
                                    
                                    space_buf(en_wbuf,sept_class) = space_buf(en_wbuf,sept_class) - 1; % remove from buffer
                                    pentry = pie{ist,sept_class};
                                    for kentry=1:K(sept_class)
                                        space_srv(en_wbuf,Ks(sept_class)+kentry) = space_srv(en_wbuf,Ks(sept_class)+kentry) + 1; % bring job in service
                                        if isSimulation
                                            % break the tie
                                            outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                            rate_k = rate;
                                            rate_k(en,:) = rate_k(en,:) * pentry(kentry);
                                            outrate = [outrate; rate_k(en,:)];
                                            outprob = [outprob; ones(size(rate(en,:),1),1)];
                                        else
                                            outspace = [outspace; space_buf(en,:), space_srv(en,:), space_var(en,:)];
                                            rate_k = rate;
                                            rate_k(en,:) = rate_k(en,:) * pentry(kentry);
                                            outrate = [outrate; rate_k(en,:)];
                                            outprob = [outprob; ones(size(rate(en,:),1),1)];
                                        end
                                        space_srv(en_wbuf,Ks(sept_class)+kentry) = space_srv(en_wbuf,Ks(sept_class)+kentry) - 1; % bring job in service
                                    end
                                otherwise
                                    line_error(mfilename,'Scheduling strategy %s is not supported.', qn.sched(ist));
                            end
                        end
                    end
                    if isSimulation
                        if size(outspace,1) > 1
                            tot_rate = sum(outrate);
                            cum_rate = cumsum(outrate) / tot_rate;
                            firing_ctr = 1 + max([0,find( rand > cum_rate' )]); % select action
                            outspace = outspace(firing_ctr,:);
                            outrate = sum(outrate);
                            outprob = outprob(firing_ctr,:);
                        end
                    end
                end
            end
        case EventType.PHASE
            outspace = [];
            outrate = [];
            outprob = [];
            [ni,nir,~,kir] = State.toMarginal(qn,ind,inspace,K,Ks,space_buf,space_srv,space_var);
            if nir(class)>0
                for k=1:(K(class)-1)
                    en = space_srv(:,Ks(class)+k) > 0;
                    if any(en)
                        for kdest=(k+1):K(class) % new phase
                            rate = 0;
                            space_srv_k = space_srv(en,:);
                            space_buf_k = space_buf(en,:);
                            space_var_k = space_var(en,:);
                            space_srv_k(:,Ks(class)+k) = space_srv_k(:,Ks(class)+k) - 1;
                            space_srv_k(:,Ks(class)+kdest) = space_srv_k(:,Ks(class)+kdest) + 1;
                            switch qn.schedid(ist)
                                case SchedStrategy.ID_EXT
                                    rate = ph{ist,class}{1}(k,kdest); % move next job forward
                                case SchedStrategy.ID_INF
                                    rate = ph{ist,class}{1}(k,kdest)*kir(:,class,k); % assume active
                                case SchedStrategy.ID_PS
                                    rate = ph{ist,class}{1}(k,kdest)*kir(:,class,k)./ni(:).*min(ni(:),S(ist)); % assume active
                                case SchedStrategy.ID_DPS
                                    if S(ist) > 1
                                        line_error(mfilename,'Multi-server DPS not supported yet');
                                    end
                                    w_i = qn.schedparam(ist,:);
                                    w_i = w_i / sum(w_i);
                                    rate = ph{ist,class}{1}(k,kdest)*kir(:,class,k)*w_i(class)./(sum(repmat(w_i,size(nir,1),1)*nir',2)); % assume active
                                case SchedStrategy.ID_GPS
                                    if S(ist) > 1
                                        line_error(mfilename,'Multi-server GPS not supported yet');
                                    end
                                    cir = min(nir,ones(size(nir)));
                                    w_i = qn.schedparam(ist,:); w_i = w_i / sum(w_i);
                                    rate = ph{ist,class}{1}(k,kdest)*kir(:,class,k)/nir(class)*w_i(class)/(w_i*cir(:)); % assume active
                                    
                                case {SchedStrategy.ID_FCFS, SchedStrategy.ID_HOL, SchedStrategy.ID_LCFS, SchedStrategy.ID_SIRO, SchedStrategy.ID_SEPT, SchedStrategy.ID_LEPT}
                                    rate = ph{ist,class}{1}(k,kdest)*kir(:,class,k); % assume active
                            end
                            % if the class cannot be served locally,
                            % then rate = NaN since mu{i,class}=NaN
                            outrate = [outrate; rate];
                            outspace = [outspace; space_buf_k, space_srv_k, space_var_k];
                            outprob = [outprob; ones(size(rate,1),1)];
                        end
                    end
                end
                if isSimulation
                    if size(outspace,1) > 1
                        tot_rate = sum(outrate);
                        cum_rate = cumsum(outrate) / tot_rate;
                        firing_ctr = 1 + max([0,find( rand > cum_rate' )]); % select action
                        outspace = outspace(firing_ctr,:);
                        outrate = sum(outrate);
                        outprob = outprob(firing_ctr,:);
                    end
                end
            end
    end
elseif qn.isstateful(ind)
    switch qn.nodetype(ind)
        case NodeType.Cache
            % job arrives in class, then reads and moves into hit or miss
            % class, then departs
            switch event
                case EventType.ARV
                    space_srv(:,class) = space_srv(:,class) + 1;
                    outspace = [space_srv, space_var]; % buf is empty
                    outrate = -1*ones(size(outspace,1)); % passive action, rate is unspecified
                case EventType.DEP
                    if space_srv(class)>0
                        space_srv(:,class) = space_srv(:,class) - 1;
                        outspace = [space_srv, space_var]; % buf is empty
                        outrate = Distrib.InfRate*ones(size(outspace,1)); % passive action, rate is unspecified
                    end
                case EventType.READ
                    n = qn.varsparam{ind}.nitems; % n items
                    m = qn.varsparam{ind}.cap; % capacity
                    ac = qn.varsparam{ind}.accost; % access cost                    
                    hitclass = qn.varsparam{ind}.hitclass;
                    missclass = qn.varsparam{ind}.missclass;
                    h = length(m);
                    rpolicy_id = qn.varsparam{ind}.rpolicy;
                    if space_srv(class)>0 && sum(space_srv)==1 % is a job of class is in
                        p = qn.varsparam{ind}.pref{class};
                        en =  space_srv(:,class) > 0;
                        space_srv_k = [];
                        space_var_k = [];
                        outrate = [];
                        if any(en)
                            for e=find(en)'
                                if isSimulation
                                    % pick one item
                                    kset = 1 + max([0,find( rand > cumsum(p) )]);
                                    % pick one entry list
                                    l = 1 + max([0,find( rand > cumsum(ac{class,kset}(1,:)) )]);
                                    outprob = ac{class,kset}(1,l) * p(kset);
                                else
                                    kset = 1:n;
                                end
                                for k=kset % request to item k
                                    space_srv_e = space_srv(e,:);
                                    space_srv_e(class) = space_srv_e(class) - 1;
                                    var = space_var(e,:);
                                    posk = find(k==var);
                                    
                                    if isempty(posk) % CACHE MISS, add to list 1
                                        space_srv_e(missclass(class)) = space_srv_e(missclass(class)) + 1;
                                        switch rpolicy_id
                                            case {ReplacementStrategy.ID_FIFO, ReplacementStrategy.ID_LRU, ReplacementStrategy.ID_SFIFO}
                                                varp = var;
                                                varp(2:m(1)) = var(1:(m(1)-1));
                                                varp(1) = k;
                                                space_srv_k = [space_srv_k; space_srv_e];
                                                space_var_k = [space_var_k; varp];
                                                if isSimulation
                                                    %% no p(k) weighting since that goes in the outprob vec
                                                    outrate(end+1,1) = Distrib.InfRate;
                                                else
                                                    for l=2
                                                        outrate(end+1,1) = ac{class,k}(1,l) * p(k) * Distrib.InfRate;
                                                    end
                                                end
                                            case ReplacementStrategy.ID_RR
                                                if isSimulation
                                                    varp = var;
                                                    r = randi(m(1),1,1);
                                                    varp(r) = k;
                                                    space_srv_k = [space_srv_k; space_srv_e];
                                                    space_var_k = [space_var_k; (varp)];
                                                    outrate(end+1,1) = Distrib.InfRate;
                                                else
                                                    for r=1:m(1) % random position in list 1
                                                        varp = var;
                                                        varp(r) = k;
                                                        space_srv_k = [space_srv_k; space_srv_e];
                                                        space_var_k = [space_var_k; (varp)];
                                                        for l=2
                                                            outrate(end+1,1) = ac{class,k}(1,l) * p(k)/m(1) * Distrib.InfRate;
                                                        end
                                                    end
                                                end
                                        end
                                    elseif posk <= sum(m(1:h-1)) % CACHE HIT in list i < h, move to list i+1
                                        space_srv_e(hitclass(class)) = space_srv_e(hitclass(class)) + 1;
                                        i = min(find(posk <= cumsum(m)));
                                        j = posk - sum(m(1:i-1));
                                        
                                        switch rpolicy_id
                                            case ReplacementStrategy.ID_FIFO
                                                if isSimulation
                                                    varp = var;
                                                    inew = i-1+probchoose(ac{class,k}(i,i:end)/sum(ac{class,k}(i,i:end)))-1; % can choose i
                                                    if inew~=i
                                                        varp(cpos(i,j)) = var(cpos(inew,m(inew)));
                                                        varp(cpos(inew,2):cpos(inew,m(inew))) = var(cpos(inew,1):cpos(inew,m(inew)-1));
                                                        varp(cpos(inew,1)) = k;
                                                    end
                                                    %varp(cpos(i,j)) = var(cpos(i+1,m(i+1)));
                                                    %varp(cpos(i+1,2):cpos(i+1,m(i+1))) = var(cpos(i+1,1):cpos(i+1,m(i+1)-1));
                                                    %varp(cpos(i+1,1)) = k;
                                                                                                        
                                                    space_srv_k = [space_srv_k; space_srv_e];
                                                    space_var_k = [space_var_k; varp];
                                                    outrate(end+1,1) = Distrib.InfRate;
                                                else
                                                    for inew = i:h
                                                        varp = var;
                                                        varp(cpos(i,j)) = var(cpos(inew,m(inew)));
                                                        varp(cpos(inew,2):cpos(inew,m(inew))) = var(cpos(inew,1):cpos(inew,m(inew)-1));
                                                        varp(cpos(inew,1)) = k;
                                                        space_srv_k = [space_srv_k; space_srv_e];
                                                        space_var_k = [space_var_k; varp];
                                                        outrate(end+1,1) = ac{class,k}(1+i,1+inew) * p(k) * Distrib.InfRate;
                                                    end
                                                end
                                            case ReplacementStrategy.ID_RR
                                                if isSimulation
                                                    inew = i-1+probchoose(ac{class,k}(i,i:end)/sum(ac{class,k}(i,i:end)))-1; % can choose i
                                                    varp = var;
                                                    r = randi(m(inew),1,1);
                                                    varp(cpos(i,j)) = var(cpos(inew,r));
                                                    varp(cpos(inew,r)) = k;
                                                    outprob = outprob * ac{class,k}(1+i,1+inew);
                                                    space_srv_k = [space_srv_k; space_srv_e];
                                                    space_var_k = [space_var_k; varp];
                                                    outrate(end+1,1) = Distrib.InfRate;
                                                else
                                                    for inew = i:h
                                                        for r=1:m(inew) % random position in new list
                                                            varp = var;
                                                            varp(cpos(i,j)) = var(cpos(inew,r));
                                                            varp(cpos(inew,r)) = k;
                                                            space_srv_k = [space_srv_k; space_srv_e];
                                                            space_var_k = [space_var_k; varp];
                                                            outrate(end+1,1) = ac{class,k}(1+i,1+inew) * p(k)/m(inew) * Distrib.InfRate;
                                                        end
                                                    end
                                                end
                                            case {ReplacementStrategy.ID_LRU, ReplacementStrategy.ID_SFIFO}
                                                if isSimulation
                                                    varp = var;
                                                    inew = i-1+ probchoose(ac{class,k}(i,i:end)/sum(ac{class,k}(i,i:end)))-1; % can choose i
                                                    varp(cpos(i,2):cpos(i,j)) = var(cpos(i,1):cpos(i,j-1));
                                                    varp(cpos(i,1)) = var(cpos(inew,m(inew)));
                                                    varp(cpos(inew,2):cpos(inew,m(inew))) = var(cpos(inew,1):cpos(inew,m(inew)-1));
                                                    varp(cpos(inew,1)) = k;
                                                    space_srv_k = [space_srv_k; space_srv_e];
                                                    space_var_k = [space_var_k; varp];
                                                    outrate(end+1,1) = Distrib.InfRate;
                                                else
                                                    for inew = i:h
                                                        varp = var;
                                                        varp(cpos(i,2):cpos(i,j)) = var(cpos(i,1):cpos(i,j-1));
                                                        varp(cpos(i,1)) = var(cpos(inew,m(inew)));
                                                        varp(cpos(inew,2):cpos(inew,m(inew))) = var(cpos(inew,1):cpos(inew,m(inew)-1));
                                                        varp(cpos(inew,1)) = k;
                                                        space_srv_k = [space_srv_k; space_srv_e];
                                                        space_var_k = [space_var_k; varp];
                                                        outrate(end+1,1) = ac{class,k}(1+i,1+inew) * p(k) * Distrib.InfRate;
                                                    end
                                                end
                                        end
                                    else % CACHE HIT in list h
                                        space_srv_e(hitclass(class)) = space_srv_e(hitclass(class)) + 1;
                                        i=h;
                                        j = posk - sum(m(1:i-1));
                                        switch rpolicy_id
                                            case ReplacementStrategy.ID_RR
                                                space_srv_k = [space_srv_k; space_srv_e];
                                                space_var_k = [space_var_k; (var)];
                                                if isSimulation
                                                    outrate(end+1,1) = Distrib.InfRate;
                                                else                                                    
                                                    outrate(end+1,1) = ac{class,k}(1+h,1+h) * p(k) * Distrib.InfRate;
                                                end
                                            case {ReplacementStrategy.ID_FIFO, ReplacementStrategy.ID_SFIFO}
                                                space_srv_k = [space_srv_k; space_srv_e];
                                                space_var_k = [space_var_k; var];
                                                if isSimulation
                                                    outrate(end+1,1) = Distrib.InfRate;
                                                else
                                                    outrate(end+1,1) = ac{class,k}(1+h,1+h) * p(k) * Distrib.InfRate;
                                                end
                                            case ReplacementStrategy.ID_LRU
                                                varp = var;
                                                varp(cpos(h,2):cpos(h,j)) = var(cpos(h,1):cpos(h,j-1));
                                                varp(cpos(h,1)) = var(cpos(h,j));
                                                space_srv_k = [space_srv_k; space_srv_e];
                                                space_var_k = [space_var_k; varp];
                                                if isSimulation
                                                    outrate(end+1,1) = Distrib.InfRate;
                                                else
                                                    outrate(end+1,1) = ac{class,k}(1+h,1+h) * p(k) * Distrib.InfRate;
                                                end
                                        end
                                    end
                                end
                            end
                            % if state is unchanged, still add with rate 0
                            outspace = [space_srv_k, space_var_k];
                        end
                    end
            end
    end % switch nodeType
end
    function pos = cpos(i,j)
        % POS = CPOS(I,J)
        
        pos = sum(m(1:i-1)) + j;
    end
end
