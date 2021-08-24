function runtime = run(self)
% TSIM = RUN()

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

Tstart=tic;
if ~exist('options','var')
    options = self.getOptions;
end

if ~isfield(options,'verbose')
    options.verbose = 0;
end

if ~isfield(options,'force')
    options.force = false;
end

if ~isfield(options,'keep')
    options.keep = false;
end

if self.enableChecks && ~self.supports(self.model)
    %    if options.verbose
    line_error(mfilename,'This model contains features not supported by the JMT solver.');
    %    end
    %    runtime = toc(T0);
    %    return
end

if ~isfield(options,'samples')
    options.samples = 1e4; % default: this is the samples / measure, not the total number of simulation events, which can be much larger.
elseif options.samples < 5e3
    if ~strcmpi(options.method,'jmva.ls')
        line_warning(mfilename,'JMT requires at least 5000 samples for each metric, the current value is %d. Starting the simulation with 5000 samples.', options.samples);
    end
    options.samples = 5e3;
end

if ~isfield(options,'verbose')
    options.verbose = 0;
end

if ~isfield(options,'keep')
    options.verbose = false;
end

if ~isfield(options,'seed')
    options.seed = randi([1,1e6]);
end
self.seed = options.seed;

if ~isfield(options,'timespan')
    options.timespan = [0,Inf];
else
    self.maxSimulatedTime = options.timespan(2);
end

if ~self.model.hasInitState
    self.model.initDefault;
end

self.maxSamples = options.samples;

switch options.method
    case {'jsim','default'}
        if true %isinf(options.timespan(2)) || ((options.timespan(2)) == (options.timespan(1)))
            self.writeJSIM;
            cmd = ['java -cp "',self.getJMTJarPath(),filesep,'JMT.jar" jmt.commandline.Jmt sim "',self.getFilePath(),'jsimg',filesep,self.getFileName(),'.jsimg" -seed ',num2str(options.seed),' --illegal-access=permit'];
            if options.verbose
                line_printf('\nJMT model: %s',[self.getFilePath(),'jsimg',filesep,self.getFileName(),'.jsimg']);                
                line_printf('\nJMT command: %s',cmd);
            end
            [~, result] = system(cmd);
            runtime = toc(Tstart);
            self.getResults;
            if ~options.keep
                delete([self.getFilePath(),'jsimg',filesep,self.getFileName(),'*']);
            end
        else
            options = self.getOptions;
            initSeed = self.options.seed;
            initTimeSpan = self.options.timespan;
            self.options.timespan(1) = self.options.timespan(2);
            if isfield(options,'timespan')  && isfinite(options.timespan(2))
                qn = self.getStruct;
                tu = [];
                for it=1:options.iter_max
                    self.options.seed = initSeed + it -1;
                    TranSysStateAggr{it} = self.sampleSysAggr();
                    if isempty(tu)
                        tu = TranSysStateAggr{it}.t;
                    else
                        % we need to limit the time series at the minimum
                        % as otherwise the predictor of the state cannot
                        % take into account constraints that exist on the
                        % state space
                        tumax = min(max(tu),max(TranSysStateAggr{it}.t));
                        tu = union(tu, TranSysStateAggr{it}.t);
                        tu = tu(tu<=tumax);
                    end
                end
                QNt = cellzeros(qn.nstations, qn.nclasses, length(tu), 2);
                UNt = cellzeros(qn.nstations, qn.nclasses, length(tu), 2);
                TNt = cellzeros(qn.nstations, qn.nclasses, length(tu), 2);
                M = qn.nstations;
                K = qn.nclasses;
                for j=1:M
                    for r=1:K
                        QNt{j,r}(:,2) = tu;
                        UNt{j,r}(:,2) = tu;
                        TNt{j,r}(:,2) = tu;
                        for it=1:options.iter_max
                            qlenAt_t = interp1(TranSysStateAggr{it}.t, TranSysStateAggr{it}.state{j}(:,r), tu,'previous');
                            avgQlenAt_t = qlenAt_t;
                            %avgQlenAt_t = cumsum(qlenAt_t .*[0;diff(tu)])./tu;
                            avgQlenAt_t(isnan(avgQlenAt_t))=0;
                            QNt{j,r}(:,1) = QNt{j,r}(:,1) + (1/options.iter_max) * avgQlenAt_t;
                        end
                        for it=1:options.iter_max
                            if isfinite(qn.nservers(j))
                                occupancyAt_t = interp1(TranSysStateAggr{it}.t, min(TranSysStateAggr{it}.state{j}(:,r),qn.nservers(j)), tu,'previous')/qn.nservers(j);
                            else % if delay we use queue-length
                                occupancyAt_t = interp1(TranSysStateAggr{it}.t, TranSysStateAggr{it}.state{j}(:,r), tu,'previous');
                            end
                            avgOccupancyAt_t = occupancyAt_t;
                            %avgOccupancyAt_t = cumsum(occupancyAt_t .*[0;diff(tu)])./tu;
                            avgOccupancyAt_t(isnan(avgOccupancyAt_t))=0;
                            UNt{j,r}(:,1) = UNt{j,r}(:,1) + (1/options.iter_max) * avgOccupancyAt_t;
                        end
                        %                         for it=1:options.iter_max
                        %                             departures = [0;diff(TranSysStateAggr{it}.state{j}(:,r))];
                        %                             departures(departures>0) = 0;
                        %                             departuresAt_t = abs(interp1(TranSysStateAggr{it}.t, cumsum(departures), tu, 'previous'));
                        %                             avgDeparturesAt_t = departuresAt_t./tu;
                        %                             avgDeparturesAt_t(isnan(avgDeparturesAt_t))=0;
                        %                             TNt{j,r}(:,1) = TNt{j,r}(:,1) + (1/options.iter_max) * avgDeparturesAt_t;
                        %                         end
                        if isfinite(qn.nservers(j))
                            TNt{j,r}(:,1) = UNt{j,r}(:,1) * qn.nservers(j) * qn.rates(j,r);
                        else
                            TNt{j,r}(:,1) = UNt{j,r}(:,1) * qn.rates(j,r);
                        end
                    end
                end
                runtime = toc(Tstart);
                RNt = [];
                CNt = [];
                XNt = [];
                self.setTranAvgResults(QNt,UNt,RNt,TNt,CNt,XNt,runtime);
                self.result.Tran.Avg.U = UNt;
                self.result.Tran.Avg.T = TNt;
                self.result.Tran.Avg.Q = QNt;
            end
            self.options.seed = initSeed;
            self.options.timespan = initTimeSpan;
            self.result.('solver') = self.getName();
            self.result.runtime = runtime;
        end
    case {'jmva','jmva.amva','jmva.mva','jmva.recal','jmva.comom','jmva.chow','jmva.bs','jmva.aql','jmva.lin','jmva.dmlin','jmva.ls',...
          'jmt.jmva','jmt.jmva.mva','jmt.jmva.amva','jmt.jmva.recal','jmt.jmva.comom','jmt.jmva.chow','jmt.jmva.bs','jmt.jmva.aql','jmt.jmva.lin','jmt.jmva.dmlin','jmt.jmva.ls'}
        fname = self.writeJMVA([self.getFilePath(),'jmva',filesep,self.getFileName(),'.jmva']);
        cmd = ['java -cp "',self.getJMTJarPath(),filesep,'JMT.jar" jmt.commandline.Jmt mva "',fname,'" -seed ',num2str(options.seed),' --illegal-access=permit'];
        if options.verbose
            line_printf('\nJMT model: %s',[self.getFilePath(),'jmva',filesep,self.getFileName(),'.jmva']);
            line_printf('\nJMT command: %s',cmd);
        end
        [~, result] = system(cmd);
        runtime = toc(Tstart);
        self.getResults;
        if ~options.keep
            delete([self.getFilePath(),'jmva',filesep,self.getFileName(),'*']);
        end
    otherwise
        line_warning(mfilename,'This solver does not support the specified method. Setting to default.');
        self.options.method  = 'default';
        runtime = run(self);
end

if options.verbose > 0
    line_printf('\nJMT analysis completed. Runtime: %f seconds.\n',runtime);
end
end
