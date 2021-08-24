classdef NetworkStruct <handle
    % Data structure representation for a Network object
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties
        cap;     % total buffer size
        chains;     % binary CxK matrix where 1 in entry (i,j) indicates that class j is in chain i.
        classcap;    % buffer size for each class
        classnames;  % name of each job class
        classprio;       % scheduling priorities in each class (optional)
        csmask; % (r,s) entry if class r can switch into class s somewhere
        %forks;      % forks table from each station
        % (MKxMK matrix with integer entries), indexed first by
        % station, then by class
        isstatedep; % state dependent routing
        isstation; % element i is true if node i is a station
        isstateful; % element i is true if node i is stateful
        isslc; % element r is true if class r self-loops at its reference station
        lst; % laplace-stieltjes transform
        mu;          % service rate in each service phase, for each job class in each station
        % (MxK cell with n_{i,k}x1 double entries)
        nchains;           % number of chains (int)
        nclasses;          % number of classes (int)
        nclosedjobs;          % total population (int)
        njobs;             % initial distribution of jobs in classes (Kx1 int)
        nnodes; % number of nodes (Mn int)
        nservers;   % number of servers per station (Mx1 int)
        nstations;  % number of stations (int)
        nstateful;  % number of stations (int)
        nvars; % number of local variables
        nodenames;   % name of each node
        nodevisits;  % visits placed by classes at the nodes
        nodetype; % server type in each node
        phases; % number of phases in each service or arrival process
        phasessz; % number of phases
        phaseshift; % shift for phases
        pie;        % probability of entry in each each service phase
        phi;         % probability of service completion in each service phase,
        % for each job class in each station
        % (MxK cell with n_{i,k}x1 double entries)
        proc;     % cell matrix of service and arrival process representations
        rates;       % service rate for each job class in each station
        refstat;    % index of the reference node for each request class (Kx1 int)
        routing;     % routing strategy type
        rt;         % routing table with class switching
        % (M*K)x(M*K) matrix with double entries), indexed first by
        % station, then by class
        rtnodes;         % routing table with class switching
        % (Mn*K)x(Mn*K) matrix with double entries), indexed first by
        % node, then by class
        rtfun; % local routing functions
        % (Mn*K)x(Mn*K) matrix with double entries), indexed first by
        % station, then by class
        sched;       % scheduling strategy in each station
        schedid;       % scheduling strategy id in each station (optional)
        schedparam;       % scheduling weights in each station and class (optional)
        sync;
        space;    % state space
        state;    % initial or current state
        scv; % squared coefficient of variation of service times (MxK)
        visits;           % visits placed by classes at the resources
        varsparam;     % parameters for local variables
    end
    
    methods
        %constructor
        function self = NetworkStruct(nodetype, nodenames, classnames, numservers, njobs, refstat, routing)
            % SELF = NETWORKSTRUCT(NODETYPE, NODENAMES, CLASSNAMES, NUMSERVERS, NJOBS, REFSTAT, ROUTING)
            self.nnodes = numel(nodenames);
            self.nclasses = length(classnames);
            self.nclosedjobs = sum(njobs(isfinite(njobs)));
            self.nservers = numservers;
            self.nodetype = -1*ones(self.nstations,1);
            self.scv = ones(self.nstations,self.nclasses);
            %self.forks = zeros(M,K);
            self.njobs = njobs(:)';
            self.refstat = refstat;
            self.space = cell(self.nstations,1);
            self.routing = routing;
            self.chains = [];
            self.lst = {};
            if exist('nvars','var') && ~isempty(nvars)
                self.nvars = nvars;
            end
            if exist('nodetype','var') && ~isempty(nodetype)
                self.nodetype = nodetype;
                self.isstation = NodeType.isStation(nodetype);
                self.nstations = sum(self.isstation);
                self.isstateful = NodeType.isStateful(nodetype);
                self.isstatedep = false(self.nnodes,3); % col 1: buffer, col 2: srv, col 3: routing
                for ind=1:self.nnodes
                    switch self.nodetype(ind)
                        case NodeType.Cache
                            self.isstatedep(ind,2) = true; % state dependent service
                            %                            self.isstatedep(ind,3) = true; % state dependent routing
                    end
                    for r=1:self.nclasses
                        switch self.routing(ind,r)
                            case {RoutingStrategy.ID_RRB, RoutingStrategy.ID_JSQ}
                                self.isstatedep(ind,3) = true; % state dependent routing
                        end
                    end
                end
                self.nstateful = sum(self.isstateful);
                self.state = cell(self.nstations,1); for i=1:self.nstateful self.state{i} = []; end
            end
            if exist('nodenames','var')
                self.nodenames = nodenames;
            else
                self.nodenames = cell(self.nnodes,1);
                for j = 1:self.nstations
                    self.nodenames{j,1} = int2str(j);
                end
            end
            if exist('classnames','var')
                self.classnames = classnames;
            else
                self.classnames = cell(self.nclasses,1);
                for j = 1:self.nclasses
                    self.classnames{j,1} = int2str(j);
                end
            end
            self.reindex();
        end
        
        function reindex(self)
            % REINDEX()
            
            for ind=1:self.nnodes
                self.nodeToStateful(ind) = self.nd2sf(ind);
                self.nodeToStation(ind) = self.nd2st(ind);
            end
            for ist=1:self.nstations
                self.stationToNode(ist) = self.st2nd(ist);
                self.stationToStateful(ist) = self.st2sf(ist);
            end
            for isf=1:self.nstateful
                self.statefulToNode(isf) = self.sf2nd(isf);
            end
        end
        
        function setChains(self, chains, visits, rt, nodes_visits)
            % SETCHAINS(CHAINS, VISITS, RT, NODES_VISITS)
            
            self.chains = logical(chains);
            self.visits = visits;
            self.rt = rt;
            self.nchains = size(chains,1);
            self.nodevisits = nodes_visits;
            self.isslc = false(self.nchains,1);
            for c=1:self.nchains
                if nnz(visits{c}) == 1
                    self.isslc(c) = true;
                end
            end
        end
        
        function setSched(self, sched, schedparam)
            % SETSCHED(SCHED, SCHEDPARAM)
            
            self.sched = sched;
            self.schedparam = schedparam;
            self.schedid = zeros(self.nstations,1);
            for i=1:self.nstations
                self.schedid(i) = SchedStrategy.toId(sched(i));
            end
        end
        
        function setLSTs(self, lst)
            % SETLAPLACETRANSFORMS(LT)
            self.lst = lst;
        end
        
        function setService(self, rates, scv)
            % SETSERVICE(RATES, SCV)
            
            self.rates = rates;
            self.scv = scv;
        end
        
        function setMAPService(self, map, phases)
            % SETMAPSERVICE(MAP, PHASES)
            
            self.proc = map;
            self.pie = cell(size(map));
            for i=1:size(map,1)
                for r=1:size(map,2)
                    if ~isempty(map{i,r})
                        self.proc{i,r} = map_normalize(map{i,r});
                        self.pie{i,r} = map_pie(map{i,r});
                    else
                        self.pie{i,r} = NaN;
                    end
                end
            end
            self.phases = phases;
            self.phasessz = max(self.phases,ones(size(self.phases)));
            self.phasessz(self.nodeToStation(self.nodetype == NodeType.Join),:)=phases(self.nodeToStation(self.nodetype == NodeType.Join),:);
            self.phaseshift = [zeros(size(phases,1),1),cumsum(self.phasessz,2)];
        end
        
        function setPHService(self, ph, phases)
            % SETPHSERVICE(PH, PHASES)
            
            self.proc = ph;
            self.pie = cell(size(ph));
            for i=1:size(ph,1)
                for r=1:size(ph,2)
                    if ~isempty(ph{i,r})
                        self.proc{i,r} = map_normalize(ph{i,r});
                        self.pie{i,r} = map_pie(ph{i,r});
                    else
                        self.pie{i,r} = NaN;
                    end
                end
            end
            self.phases = phases;
            self.phasessz = max(self.phases,ones(size(self.phases)));
            self.phasessz(self.nodeToStation(self.nodetype == NodeType.Join),:)=phases(self.nodeToStation(self.nodetype == NodeType.Join),:);
            self.phaseshift = [zeros(size(phases,1),1),cumsum(self.phasessz,2)];
        end
        
        function setCoxService(self, mu, phi, phases)
            % SETCOXSERVICE(MU, PHI, PHASES)
            
            self.mu = mu;
            self.phi = phi;
            self.phases = phases;
            self.phasessz = max(self.phases,ones(size(self.phases)));
            self.phaseshift = [zeros(size(phases,1),1),cumsum(self.phasessz,2)];
        end
        
        function setCapacity(self, cap, classcap)
            % SETCAPACITY(CAP, CLASSCAP)
            
            self.cap = cap;
            self.classcap = classcap;
        end
        
        function setPrio(self, prio)
            % SETPRIO(PRIO)
            
            self.classprio = prio;
        end
        
        function setSync(self, sync)
            % SETSYNC(SYNC)
            
            self.sync = sync;
        end
        
        function setLocalVars(self, nvars, varsparam)
            % SETLOCALVARS(NVARS, VARSPARAM)
            
            self.nvars = nvars;
            self.varsparam = varsparam;
        end
        
        function setRoutingFunction(self, rtfun, csmask)
            %(RTFUN, CSMASK)
            
            self.rtfun = rtfun;
            self.csmask = csmask;
        end
    end
    
    properties (Access = public)
        nodeToStateful;
        nodeToStation;
        stationToNode;
        stationToStateful;
        statefulToNode;
    end
    
    methods % index conversion
        function stat_idx = nd2st(self, node_idx)
            % STAT_IDX = ND2ST(NODE_IDX)
            
            if self.isstation(node_idx)
                stat_idx = at(cumsum(self.isstation),node_idx);
            else
                stat_idx = NaN;
            end
        end
        
        function node_idx = st2nd(self,stat_idx)
            % NODE_IDX = ST2ND(SELF,STAT_IDX)
            
            v = cumsum(self.isstation) == stat_idx;
            if any(v)
                node_idx =  find(v, 1);
            else
                node_idx = NaN;
            end
        end
        
        function sful_idx = st2sf(self,stat_idx)
            % SFUL_IDX = ST2SF(SELF,STAT_IDX)
            
            sful_idx = nd2sf(self,st2nd(self,stat_idx));
        end
        
        function sful_idx = nd2sf(self, node_idx)
            % SFUL_IDX = ND2SF(NODE_IDX)
            
            if self.isstateful(node_idx)
                sful_idx = at(cumsum(self.isstateful),node_idx);
            else
                sful_idx = NaN;
            end
        end
        
        function node_idx = sf2nd(self,stat_idx)
            % NODE_IDX = SF2ND(SELF,STAT_IDX)
            
            v = cumsum(self.isstateful) == stat_idx;
            if any(v)
                node_idx =  find(v, 1);
            else
                node_idx = NaN;
            end
        end
    end
    
    methods (Access = public)
        function bool = isclosed(self)
            bool = all(isfinite(self.njobs));
        end
        
        function bool = isopen(self)
            bool = all(isinf(self.njobs));
        end
        
        function bool = ismixed(self)
            bool = true;
            if isopen(self) || isclosed(self)
                bool = false;
            end
        end
        
        function bool = isqsys(self)
            % bool = isqsys(self)
            % True if the model defined an open queueing system.
            bool = false;
            if qn.nstations == 2 && sum(self.nodetype == NodeType.Sink)==1
                bool = true;
            end
        end
        
        function bool = hasdelay(self)
            bool = any(self.nodetype == NodeType.Delay);
        end
        
        function newObj = copy(obj)
            % NEWOBJ = COPY(OBJ)
            
            newObj = NetworkStruct(obj.nodetype, obj.nodenames, obj.classnames, obj.nservers, obj.njobs, obj.refstat, obj.routing);
            newObj.cap = obj.cap;     % total buffer size
            newObj.chains = obj.chains;     % binary CxK matrix where 1 in entry (i,j) indicates that class j is in chain i.
            newObj.classcap = obj.classcap;    % buffer size for each class
            newObj.classnames = obj.classnames;  % name of each job class
            newObj.classprio = obj.classprio;       % scheduling priorities in each class (optional)
            newObj.csmask = obj.csmask; % (r,s) entry if class r can switch into class s somewhere
            newObj.isstatedep = obj.isstatedep; % state dependent routing
            newObj.isstation = obj.isstation; % element i is true if node i is a station
            newObj.isstateful = obj.isstateful; % element i is true if node i is stateful
            newObj.isslc = obj.isslc; % element r is true if class r self-loops at its reference station
            newObj.lst = obj.lst;          % function handle to Laplace-Stieltjes transform for each service or arrival distribution
            newObj.mu = obj.mu;          % service rate in each service phase, for each job class in each station
            newObj.nchains = obj.nchains;           % number of chains (int)
            newObj.nclasses = obj.nclasses;          % number of classes (int)
            newObj.nclosedjobs = obj.nclosedjobs;          % total population (int)
            newObj.njobs = obj.njobs;             % initial distribution of jobs in classes (Kx1 int)
            newObj.nnodes = obj.nnodes; % number of nodes (Mn int)
            newObj.nservers = obj.nservers;   % number of servers per station (Mx1 int)
            newObj.nstations = obj.nstations;  % number of stations (int)
            newObj.nstateful = obj.nstateful;  % number of stations (int)
            newObj.nvars = obj.nvars; % number of local variables
            newObj.nodenames = obj.nodenames;   % name of each node
            newObj.nodevisits = obj.nodevisits;   % name of each node
            newObj.nodetype = obj.nodetype; % server type in each node
            newObj.phases = obj.phases; % number of phases in each service or arrival process
            newObj.phasessz = obj.phasessz; % number of phases in each service or arrival process
            newObj.phaseshift = obj.phaseshift; % number of phases in each service or arrival process
            newObj.phi = obj.phi;         % probability of service completion in each service phase,
            newObj.pie = obj.pie;         % probability of service entry in each service phase,
            newObj.proc = obj.proc;         % probability of service completion in each service phase,
            newObj.rates = obj.rates;       % service rate for each job class in each station
            newObj.refstat = obj.refstat;    % index of the reference node for each request class (Kx1 int)
            newObj.routing = obj.routing;     % routing strategy type
            newObj.rt = obj.rt;         % routing table with class switching
            newObj.rtnodes = obj.rtnodes;         % routing table with class switching
            newObj.rtfun = obj.rtfun; % local routing functions
            newObj.sched = obj.sched;       % scheduling strategy in each station
            newObj.schedid = obj.schedid;       % scheduling strategy id in each station (optional)
            newObj.schedparam = obj.schedparam;       % scheduling weights in each station and class (optional)
            newObj.sync = obj.sync;
            newObj.space = obj.space;    % state space
            newObj.state = obj.state;    % initial or current state
            newObj.scv = obj.scv; % squared coefficient of variation of service times (MxK)
            newObj.visits = obj.visits;           % visits placed by classes at the resources
            newObj.varsparam = obj.varsparam;     % parameters for local variables
        end
    end % getIndex
    
%     methods (Static, Access = public)
%         function newObj = fromJSON(qn_json)
%             % NEWOBJ = FROMJSON(QN_JSON)
%             qn = jsondecode(qn_json);
%             if ismatrix(qn.visits)
%                 qn.visits = {qn.visits'};
%             end
%             newObj = NetworkStruct(qn.nodetype, qn.nodenames, qn.classnames, qn.nservers, qn.njobs, qn.refstat, qn.routing);
%             newObj.cap = qn.cap;     % total buffer size
%             newObj.chains = qn.chains;     % binary CxK matrix where 1 in entry (i,j) indicates that class j is in chain i.
%             newObj.classcap = qn.classcap;    % buffer size for each class
%             newObj.classnames = qn.classnames;  % name of each job class
%             newObj.classprio = qn.classprio;       % scheduling priorities in each class (optional)
%             newObj.csmask = qn.csmask; % (r,s) entry if class r can switch into class s somewhere
%             newObj.isstatedep = qn.isstatedep; % state dependent routing
%             newObj.isstation = qn.isstation; % element i is true if node i is a station
%             newObj.isstateful = qn.isstateful; % element i is true if node i is stateful
%             newObj.isslc = qn.isslc; % element r is true if class r self-loops at its reference station
%             newObj.lst = qn.lst;          % function handle to Laplace-Stieltjes transform for each service or arrival distribution
%             newObj.mu = qn.mu;          % service rate in each service phase, for each job class in each station
%             newObj.nchains = qn.nchains;           % number of chains (int)
%             newObj.nclasses = qn.nclasses;          % number of classes (int)
%             newObj.nclosedjobs = qn.nclosedjobs;          % total population (int)
%             newObj.njobs = qn.njobs;             % initial distribution of jobs in classes (Kx1 int)
%             newObj.nnodes = qn.nnodes; % number of nodes (Mn int)
%             newObj.nservers = qn.nservers;   % number of servers per station (Mx1 int)
%             newObj.nstations = qn.nstations;  % number of stations (int)
%             newObj.nstateful = qn.nstateful;  % number of stations (int)
%             newObj.nvars = qn.nvars; % number of local variables
%             newObj.nodenames = qn.nodenames;   % name of each node
%             newObj.nodetype = qn.nodetype; % server type in each node
%             newObj.phases = qn.phases; % number of phases in each service or arrival process
%             newObj.phasessz = qn.phasessz; % number of phases in each service or arrival process
%             newObj.phaseshift = qn.phaseshift; % number of phases in each service or arrival process
%             newObj.phi = qn.phi;         % probability of service completion in each service phase,
%             newObj.pie = qn.pie;         % probability of service entry in each service phase,
%             newObj.proc = qn.proc;         % probability of service completion in each service phase,
%             newObj.rates = qn.rates;       % service rate for each job class in each station
%             newObj.refstat = qn.refstat;    % index of the reference node for each request class (Kx1 int)
%             newObj.routing = qn.routing;     % routing strategy type
%             newObj.rt = qn.rt;         % routing table with class switching
%             newObj.rtnodes = qn.rtnodes;         % routing table with class switching
%             newObj.rtfun = qn.rtfun; % local routing functions
%             newObj.sched = qn.sched;       % scheduling strategy in each station
%             newObj.schedid = qn.schedid;       % scheduling strategy id in each station (optional)
%             newObj.schedparam = qn.schedparam;       % scheduling weights in each station and class (optional)
%             newObj.sync = qn.sync;
%             newObj.space = qn.space;    % state space
%             newObj.state = qn.state;    % initial or current state
%             newObj.scv = qn.scv; % squared coefficient of variation of service times (MxK)
%             newObj.visits = qn.visits;           % visits placed by classes at the resources
%             newObj.varsparam = qn.varsparam;     % parameters for local variables
%         end
%     end
end
