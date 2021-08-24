classdef SolverCTMC < NetworkSolver
    % A solver based on continuous-time Markov chain (CTMC) formalism.
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    methods
        function self = SolverCTMC(model,varargin)
            % SELF = SOLVERCTMC(MODEL,VARARGIN)            
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
        end
        
        runtime = runAnalysis(self, options, config)
        Pnir = getProb(self, node, state)
        Pn = getProbSys(self)
        Pnir = getProbAggr(self, ist)
        
        Pn = getProbSysAggr(self)
        [Pi_t, SSsysa] = getTranProbSysAggr(self)   
        [Pi_t, SSnode_a] = getTranProbAggr(self, node)
        [Pi_t, SSsys] = getTranProbSys(self)     
        [Pi_t, SSnode] = getTranProb(self, node)
        %RD = getCdfRespT(self, R)
        
        function [stateSpace,nodeStateSpace] = getStateSpace(self)
            % [STATESPACE, MARGSTATESPACE] = GETSTATESPACE()
            
            options = self.getOptions;
            if options.force
                self.runAnalysis;
            end
            if isempty(self.result) || ~isfield(self.result,'space')
                line_warning(mfilename,'The model solution is not available yet or has not been cached. Either solve it or use the ''force'' option to require this is done automatically, e.g., SolverCTMC(model,''force'',true).getStateSpace()');
                stateSpace = [];
                nodeStateSpace = [];
            else
                stateSpace = self.result.space;
                shift = 1;
                for i=1:length(self.result.nodeSpace)                    
                    nodeStateSpace{i} = self.result.space(:,shift:(shift+size(self.result.nodeSpace{i},2)-1));
                    shift = shift + size(self.result.nodeSpace{i},2);
                end
            end
        end
        
        function stateSpaceAggr = getStateSpaceAggr(self)
            % STATESPACEAGGR = GETSTATESPACEAGGR()
            
            options = self.getOptions;
            if options.force
                self.run;
            end
            if isempty(self.result) || ~isfield(self.result,'spaceAggr')
                line_warning(mfilename,'The model has not been cached. Either solve it or use the ''force'' option to require this is done automatically, e.g., SolverCTMC(model,''force'',true).getStateSpaceAggr()');
                stateSpaceAggr = [];
            else
                stateSpaceAggr = self.result.spaceAggr;
            end
        end
        
        function [infGen, eventFilt, ev] = getSymbolicGenerator(self)
            if ~isdeployed
                [~, F] = getGenerator(self);
                infGen = sym(zeros(size(F{1})));
                eventFilt = cell(1, length(F));
                for e = 1:length(F)
                    F{e} = full(F{e});
                    minF = min(min(F{e}(F{e}>0)));
                    if ~isempty(minF)
                        F{e} = F{e} / minF;
                        eventFilt{e} = F{e} * sym(['x',num2str(e)],'real');
                        infGen = infGen + eventFilt{e};
                    end
                end
                infGen = ctmc_makeinfgen(infGen);
                ev = self.model.getStruct.sync;
            else
                infGen = [];
                eventFilt = [];
                ev = [];
            end
        end
        
        function [infGen, eventFilt, ev] = getInfGen(self)
            [infGen, eventFilt, ev] = getGenerator(self);
        end
        
        function [infGen, eventFilt, ev] = getGenerator(self, force)
            % [INFGEN, EVENTFILT] = GETGENERATOR()
                        
            % [infGen, eventFilt] = getGenerator(self)
            % returns the infinitesimal generator of the CTMC and the
            % associated filtration for each event
            
            if nargin==1
                force=false;
            end
            options = self.getOptions;
            if options.force || force
                self.runAnalysis;
            end
            if isempty(self.result) || ~isfield(self.result,'infGen')
                line_warning(mfilename,'The model has not been cached. Either solve it, set options.force=true or call getGenerator(true).');
                infGen = [];
                eventFilt = [];
            else
                infGen = self.result.infGen;
                eventFilt = self.result.eventFilt;
            end
            ev = self.model.getStruct.sync;
        end
        
        tstate = sampleSys(self, numevents);
        sampleAggr = sampleAggr(self, node, numSamples);
        
    end
    
    methods (Static)
        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()
            
            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Source','Sink',...
                'ClassSwitch','DelayStation','Queue',...
                'APH','Coxian','Erlang','Exponential','HyperExp',...
                'StatelessClassSwitcher','InfiniteServer','SharedServer','Buffer','Dispatcher',...
                'Cache','CacheClassSwitcher', ...
                'Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS',...
                'SchedStrategy_DPS','SchedStrategy_GPS',...
                'SchedStrategy_RAND','SchedStrategy_SEPT',...
                'SchedStrategy_LEPT','SchedStrategy_FCFS',...
                'SchedStrategy_HOL','SchedStrategy_LCFS',...
                'RoutingStrategy_RRB',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'ClosedClass','OpenClass','Replayer'});
        end
        
        function [bool, featSupported, featUsed] = supports(model)
            % [BOOL, FEATSUPPORTED, FEATUSED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverCTMC.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end
        
        function checkOptions(options)
            % CHECKOPTIONS(OPTIONS)
            
            % do nothing
        end
        
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = lineDefaults('CTMC');
        end
        
        function printInfGen(Q,SS)
            % PRINTINFGEN(Q,SS)
            
            SS=full(SS);
            Q=full(Q);
            for s=1:size(SS,1)
                for sp=1:size(SS,1)
                    if Q(s,sp)>0
                        line_printf('\n%s->%s: %f',mat2str(SS(s,:)),mat2str(SS(sp,:)),double(Q(s,sp)));
                    end
                end
            end
        end
        
        function printEventFilt(sync,D,SS,myevents)
            % PRINTEVENTFILT(SYNC,D,SS,MYEVENTS)
            
            if ~exist('events','var')
                myevents = 1:length(sync);
            end
            SS=full(SS);
            for e=myevents
                D{e}=full(D{e});
                for s=1:size(SS,1)
                    for sp=1:size(SS,1)
                        if D{e}(s,sp)>0
                            line_printf('\n%s-- %d: (%d,%d) => (%d,%d) -->%s: %f',mat2str(SS(s,:)),e,sync{e}.active{1}.node,sync{e}.active{1}.class,sync{e}.passive{1}.node,sync{e}.passive{1}.class,mat2str(SS(sp,:)),double(D{e}(s,sp)));
                        end
                    end
                end
            end
        end
    end
end
