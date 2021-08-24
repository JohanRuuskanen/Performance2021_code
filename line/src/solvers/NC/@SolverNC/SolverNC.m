classdef SolverNC < NetworkSolver
    % A solver based on normalizing constant methods.
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    methods
        function self = SolverNC(model,varargin)
            % SELF = SOLVERNC(MODEL,VARARGIN)
            
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            persistent isNCLibLoaded;
            if isempty(isNCLibLoaded)
                javaaddpath(which('pfqn_nclib.jar'));
                isNCLibLoaded = true;
            end
        end
        
        runtime = runAnalysis(self, options, config)
        Pnir = getProb(self, node, state)
        Pnir = getProbAggr(self, node, state_a)
        Pn   = getProbSys(self)        
        Pn   = getProbSysAggr(self)
        RD = getCdfRespT(self, R);

        [lNormConst] = getProbNormConstAggr(self)
    end
    
    methods (Static)
        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()
            
            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink','Source',...
                'ClassSwitch','DelayStation','Queue',...
                'APH','Coxian','Erlang','Exponential','HyperExp',...
                'StatelessClassSwitcher','InfiniteServer',...
                'SharedServer','Buffer','Dispatcher',...
                'Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS','SchedStrategy_RAND',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'SchedStrategy_FCFS','ClosedClass',...
                'Cache','CacheClassSwitcher'});
                %'OpenClass',...
        end
        
        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverNC.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end

        function checkOptions(options)
            % CHECKOPTIONS(OPTIONS)            
            solverName = mfilename;
            if isfield(options,'timespan') && isfinite(options.timespan(2))
                line_error(mfilename,'Finite timespan not supported in %s',solverName);
            end
        end
        
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = lineDefaults('NC');            
        end
    end
end
