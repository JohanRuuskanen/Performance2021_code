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

switch options.method
    case {'default','stateindep','statedep'}
        % do nothing
    otherwise
    line_warning(mfilename,'This solver does not support the specified method. Setting to default.');
    options.method  = 'default';
end

if isinf(options.timespan(1))
    if options.verbose  == 2
        line_warning(mfilename,'%s requires options.timespan(1) to be finite. Setting it to 0.',mfilename);
    end
    options.timespan(1) = 0;
end

if options.timespan(1) == options.timespan(2)
    line_warning(mfilename,'%s: timespan is a single point, unsupported. Setting options.timespace(1) to 0.\n',mfilename);
    options.timespan(1) = 0;
end

if self.enableChecks && ~self.supports(self.model)
    %                if options.verbose
    %line_warning(mfilename,'This model contains features not supported by the solver.'); ME = MException('Line:FeatureNotSupportedBySolver', 'This model contains features not supported by the solver.'); throw(ME);
    %                end
    %runtime = toc(T0);
    %return
end

qn = self.model.getStruct().copy; % this gets modified later on so pass by copy
M = self.model.getNumberOfStations;
K = self.model.getNumberOfClasses;

%%
RT = 0;
lastSol= [];
Q = zeros(M,K); R = zeros(M,K); T = zeros(M,K);
U = zeros(M,K); C = zeros(1,K); X = zeros(1,K);
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
        [Qfull, Ufull, Rfull, Tfull, Cfull, Xfull, t, Qfull_t, Ufull_t, Tfull_t, lastSol] = solver_fluid_analysis(qn, options);
        [t,uniqueIdx] = unique(t);
        if isempty(lastSol) % if solution fails
            Q = NaN*ones(M,K); R = NaN*ones(M,K);
            T = NaN*ones(M,K); U = NaN*ones(M,K);
            C = NaN*ones(1,K); X = NaN*ones(1,K);
            Qt = cell(M,K); Ut = cell(M,K); Tt = cell(M,K);
            for ist=1:M
                for r=1:K
                    Qt{ist,r} = [NaN,NaN];
                    Ut{ist,r} = [NaN,NaN];
                    Tt{ist,r} = [NaN,NaN];
                end
            end
        else
            if isempty(self.result) && ~exist('Qt','var')
                Q = Qfull*s0prior_val;
                R = Rfull*s0prior_val;
                T = Tfull*s0prior_val;
                U = Ufull*s0prior_val;
                C = Cfull*s0prior_val;
                X = Xfull*s0prior_val;
                Qt = cell(M,K);
                Ut = cell(M,K);
                Tt = cell(M,K);
                for ist=1:M
                    for r=1:K
                        Qfull_t{ist,r} = Qfull_t{ist,r}(uniqueIdx);
                        Ufull_t{ist,r} = Ufull_t{ist,r}(uniqueIdx);
                        Tfull_t{ist,r} = Tfull_t{ist,r}(uniqueIdx);
                        Qt{ist,r} = [Qfull_t{ist,r} * s0prior_val,t];
                        Ut{ist,r} = [Ufull_t{ist,r} * s0prior_val,t];
                        Tt{ist,r} = [Tfull_t{ist,r} * s0prior_val,t];
                    end
                end
            else
                Q = Q + Qfull*s0prior_val;
                R = R + Rfull*s0prior_val;
                T = T + Tfull*s0prior_val;
                U = U + Ufull*s0prior_val;
                C = C + Cfull*s0prior_val;
                X = X + Xfull*s0prior_val;
                for ist=1:M
                    for r=1:K
                        [t,uniqueIdx] = unique(t);
                        Qfull_t{ist,r} = Qfull_t{ist,r}(uniqueIdx);
                        Ufull_t{ist,r} = Ufull_t{ist,r}(uniqueIdx);
                        %                                  Tfull_t{i,r} = Tfull_t{i,r}(uniqueIdx);
                        
                        tunion = union(Qt{ist,r}(:,2), t);
                        dataOld = interp1(Qt{ist,r}(:,2),Qt{ist,r}(:,1),tunion);
                        dataNew = interp1(t,Qfull_t{ist,r},tunion);
                        Qt{ist,r} = [dataOld + s0prior_val * dataNew, tunion];
                        
                        dataOld = interp1(Ut{ist,r}(:,2),Ut{ist,r}(:,1),tunion);
                        dataNew = interp1(t,Ufull_t{ist,r},tunion);
                        Ut{ist,r} = [dataOld + s0prior_val * dataNew, tunion];
                        
                        %                                 dataOld = interp1(Tt{i,r}(:,2),Tt{i,r}(:,1),tunion);
                        %                                 dataNew = interp1(t,Tfull_t{i,r},tunion);
                        %                                 Tt{i,r} = [dataOld + s0prior_val * dataNew, tunion];
                    end
                end
            end
        end
    end
    s0_id=pprod(s0_id,s0_sz-1); % update initial state
end
runtime = toc(T0);
self.result.solverSpecific = lastSol;
self.setAvgResults(Q,U,R,T,C,X,runtime);
Rt={}; Xt={}; Ct={};
self.setTranAvgResults(Qt,Ut,Rt,Tt,Ct,Xt,runtime);
end