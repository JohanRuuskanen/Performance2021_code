classdef Library < NetworkSolver
    % Library of methods to solve specific models.
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    methods
        function self = Library(model,varargin)
            % SELF = NETWORKSOLVERLIBRARY(MODEL,VARARGIN)
            
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            if strcmp(self.getOptions.method,'default')
                line_error(mfilename,'Line:UnsupportedMethod','This solver does not have a default solution method. Used the method option to choose a solution technique.');
            end
        end
        
        function runtime = run(self)
            % RUNTIME = RUN()
            % Run the solver
                        
            T0=tic;
            options = self.getOptions;
            
            if self.enableChecks && ~self.supports(self.model)
                %                if options.verbose
               %line_warning(mfilename,'This model contains features not supported by the solver.'); 
ME = MException('Line:FeatureNotSupportedBySolver', 'This model contains features not supported by the solver.'); 
throw(ME);
                %                end
                %                runtime = toc(T0);
                %                return
            end
            
            Solver.resetRandomGeneratorSeed(options.seed);
            
            [qn] = self.model.getStruct();
            
            [Q,U,R,T,C,X] = solver_lib_analysis(qn, options);
            
            runtime=toc(T0);
            self.setAvgResults(Q,U,R,T,C,X,runtime);
        end
    end
    
    methods (Static)
        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()
            
            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink','Source','Queue',...
                'Cox2','Erlang','Exponential','HyperExp',...
                'Pareto','Uniform','Det', ...
                'Buffer','Server','JobSink','RandomSource','ServiceTunnel',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'SchedStrategy_HOL','SchedStrategy_FCFS','OpenClass','Replayer'});
        end
        
        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = Library.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end
    end
    
    methods (Static)
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            
            options = Solver.defaultOptions();
            options.timespan = [Inf,Inf];
        end
    end
end
