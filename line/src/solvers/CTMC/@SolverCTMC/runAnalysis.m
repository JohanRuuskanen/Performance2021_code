function runtime = runAnalysis(self, options, config)
% RUNTIME = RUN()
% Run the solver

T0=tic;

if nargin<2
    options = self.getOptions;
end
if nargin<3
    config = [];
end

if ~isinf(options.timespan(1)) && (options.timespan(1) == options.timespan(2))
    line_warning(mfilename,'%s: timespan is a single point, spacing by options.tol (%e).\n',mfilename, options.tol);
    options.timespan(2) = options.timespan(1) + options.tol;
end

Solver.resetRandomGeneratorSeed(options.seed);

if self.enableChecks && ~self.supports(self.model)
    %                if options.verbose
   %line_warning(mfilename,'This model contains features not supported by the solver.'); 
ME = MException('Line:FeatureNotSupportedBySolver', 'This model contains features not supported by the solver.'); 
throw(ME);
    %                end
    %                runtime = toc(T0);
    %                return
end

qn = self.getStruct();

if any(isinf(qn.njobs))
    if isinf(options.cutoff)
        line_error(mfilename,'Line:NoCutoff','The model has open chains, it is mandatory to specify a finite cutoff value, e.g., SolverCTMC(model,''cutoff'',1).');
    end
end

M = qn.nstations;
K = qn.nclasses;
NK = qn.njobs;
sizeEstimator = 0;
for k=1:K
    sizeEstimator = sizeEstimator + gammaln(1+NK(k)+M-1) - gammaln(1+M-1) - gammaln(1+NK(k)); % worst-case estimate of the state space
end
if sizeEstimator > 6
    if ~isfield(options,'force') || options.force == false
%        line_error(mfilename,'Line:ModelTooLargeToSolve','CTMC size may be too large to solve. Stopping SolverCTMC. Set options.force=true to bypass this control.\n');
        ME = MException('Line:ModelTooLargeToSolve', 'CTMC size may be too large to solve. Stopping SolverCTMC. Set options.force=true to bypass this control.\n');
        throw(ME);
        return
    end
end

% we compute all metrics anyway because CTMC has essentially
% the same cost
if isinf(options.timespan(1))
    [s0, s0prior] = self.model.getState;
    for ind=1:qn.nnodes
        if qn.isstateful(ind)
            isf = qn.nodeToStateful(ind);
            qn.state{isf} = s0{isf}(maxpos(s0prior{1}),:); % pick one particular initial state
        end
    end
    [QN,UN,RN,TN,CN,XN,Q,SS,SSq,Dfilt,~,~,qn] = solver_ctmc_analysis(qn, options);
    % update initial state if this has been corrected by the state space
    % generator
    for isf=1:qn.nstateful
        ind = qn.statefulToNode(isf);
        self.model.nodes{ind}.setState(qn.state{isf});
        switch class(self.model.nodes{qn.statefulToNode(isf)})
            case 'Cache'
                self.model.nodes{qn.statefulToNode(isf)}.server.actualHitProb = qn.varsparam{ind}.actualhitprob;
                self.model.nodes{qn.statefulToNode(isf)}.server.actualMissProb = qn.varsparam{ind}.actualmissprob;
                self.model.refreshChains(true);                
        end
    end
    %qn.space = SS;
    self.result.infGen = Q;
    self.result.space = SS;
    self.result.spaceAggr = SSq;
    self.result.nodeSpace = qn.space;
    self.result.eventFilt = Dfilt;
    runtime = toc(T0);
    qn.space = {};
    self.setAvgResults(QN,UN,RN,TN,CN,XN,runtime,options.method);
else
    lastSol= [];
    [s0, s0prior] = self.model.getState;
    s0_sz = cellfun(@(x) size(x,1), s0)';
    s0_id = pprod(s0_sz-1);
    while s0_id>=0 % for all possible initial states
        s0prior_val = 1;
        for ind=1:qn.nnodes
            if qn.isstateful(ind)
                isf = qn.nodeToStateful(ind);
                s0prior_val = s0prior_val * s0prior{isf}(1+s0_id(isf)); % update prior
                qn.state{isf} = s0{isf}(1+s0_id(isf),:); % assign initial state to network
            end
        end
        if s0prior_val > 0
            [t,pit,QNt,UNt,~,TNt,~,~,Q,SS,SSq,Dfilt,runtime_t] = solver_ctmc_transient_analysis(qn, options);
            self.result.space = SS;
            self.result.spaceAggr = SSq;
            self.result.infGen = Q;
            self.result.eventFilt = Dfilt;
            %qn.space = SS;
            setTranProb(self,t,pit,SS,runtime_t);
            if isempty(self.result) || ~isfield(self.result,'Tran') || ~isfield(self.result.Tran,'Avg') || ~isfield(self.result.Tran.Avg,'Q')
                self.result.Tran.Avg.Q = cell(M,K);
                self.result.Tran.Avg.U = cell(M,K);
                self.result.Tran.Avg.T = cell(M,K);
                for i=1:M
                    for r=1:K
                        self.result.Tran.Avg.Q{i,r} = [QNt{i,r} * s0prior_val,t];
                        self.result.Tran.Avg.U{i,r} = [UNt{i,r} * s0prior_val,t];
                        self.result.Tran.Avg.T{i,r} = [TNt{i,r} * s0prior_val,t];
                    end
                end
            else
                for i=1:M
                    for r=1:K
                        tunion = union(self.result.Tran.Avg.Q{i,r}(:,2), t);
                        dataOld = interp1(self.result.Tran.Avg.Q{i,r}(:,2),self.result.Tran.Avg.Q{i,r}(:,1),tunion);
                        dataNew = interp1(t,QNt{i,r},tunion);
                        self.result.Tran.Avg.Q{i,r} = [dataOld+s0prior_val*dataNew,tunion];
                        dataOld = interp1(self.result.Tran.Avg.U{i,r}(:,2),self.result.Tran.Avg.U{i,r}(:,1),tunion);
                        dataNew = interp1(t,UNt{i,r},tunion);
                        self.result.Tran.Avg.U{i,r} = [dataOld+s0prior_val*dataNew,tunion];
                        
                        dataOld = interp1(self.result.Tran.Avg.T{i,r}(:,2),self.result.Tran.Avg.T{i,r}(:,1),tunion);
                        dataNew = interp1(t,TNt{i,r},tunion);
                        self.result.Tran.Avg.T{i,r} = [dataOld+s0prior_val*dataNew,tunion];
                    end
                end
            end
        end
        s0_id=pprod(s0_id,s0_sz-1); % update initial state
    end
    runtime = toc(T0);
    qn.space = {};
    self.result.('solver') = self.getName();
    self.result.runtime = runtime;
    self.result.solverSpecific = lastSol;
end
end