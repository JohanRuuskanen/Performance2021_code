classdef Solver < handle
    % Abstract class for model solution algorithms and tools
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        VERBOSE_SILENT = 0;
        VERBOSE_STD = 1;
        VERBOSE_DEBUG = 2;
    end

    properties (Hidden, Access = public)
        enableChecks;
    end
        
    properties (Access = public)
        options; % Data structure with solver options
        name; % Solver name
        model; % Model to be solved
        result; % last result
    end
    
    methods (Hidden)
        %Constructor
        function self = Solver(model, name, options)
            % SELF = SOLVER(MODEL, NAME, OPTIONS)
            if ~exist('options','var')
                options = self.defaultOptions();
            end
            self.model = model;
            self.name = name;
            self.options = options;
            self.enableChecks = true;
        end
    end
    
    methods %(Abstract) % implemented with errors for Octave compatibility
        
        function model = getModel(self)
            model = self.model;
        end
        
        function bool = supports(self,model)
            % BOOL = SUPPORTS(SELF,MODEL)
            % True if the input model is supported by the solver
            line_error(mfilename,'Line:AbstractMethodCall','An abstract method was called. The function needs to be overridden by a subclass.');            
        end
        
        function runtime = runAnalysis(self, options, config) % generic method to run the solver
            % RUNTIME = RUN()
            % Run the solver % GENERIC METHOD TO RUN THE SOLVER
            % Solve the model
            line_error(mfilename,'Line:AbstractMethodCall','An abstract method was called. The function needs to be overridden by a subclass.');                        
        end
        
        function self = setChecks(self, bool)
            self.enableChecks = bool;
        end
    end
    
    methods
        function out = getName(self)
            % OUT = GETNAME()
            % Get solver name
            out = self.name;
        end
        
        function options = getDefaultOptions(self)
            % OPTIONS = GETDEFAULTOPTIONS()
            % Get option data structure with default values
            options = Solver.defaultOptions;
        end
        
    end
    
    methods
        
        function checkOptions(self, ~)
            % CHECKOPTIONS(~)
            % Check if input option data structure is valid for the given model
            if strcmp(self.options.method,'exact')
                if ~self.model.hasProductFormSolution
                    line_error(mfilename,'The model does not have a product-form solution, hence exact solution is not possible.');
                end
            end
        end
        
        function results = getResults(self)
            % RESULTS = GETRESULTS()
            % Return results data structure
            results = self.result;
        end
        
        function bool = hasResults(self)
            % BOOL = HASRESULTS()
            % Check if the model has been solved
            bool = ~isempty(self.result);
        end
        
        function options = getOptions(self)
            % OPTIONS = GETOPTIONS()
            % Return options data structure
            options = self.options;
        end
        
		function reset(self)
            % RESET()
            % Dispose previously stored results
			self.resetResults();
		end
		
        function resetResults(self)
            % RESET()
            % Dispose previously stored results
            self.result = [];
        end
        
        function self = setOptions(self, options)
            % SELF = SETOPTIONS(OPTIONS)
            % Set a new options data structure
            defaultOptions = self.defaultOptions;
            optList = Solver.listValidOptions();
            for l=1:length(optList)
                if ~isfield(options,optList{l})
                    options.(optList{l}) = defaultOptions.(optList{l});
                end
            end
            self.options = options;
        end
    end
    
    methods (Static)
        
        function resetRandomGeneratorSeed(seed)
            % RESETRANDOMGENERATORSEED(SEED)
            % Assign a new seed to the random number generator
            if ~isoctave
                warning('off','MATLAB:RandStream:ActivatingLegacyGenerators');
                warning('off','MATLAB:RandStream:ReadingInactiveLegacyGeneratorState');
            end
            rand('seed',seed);
            %rng(seed,'twister');            
        end
        
        function bool = isAvailable()
            % BOOL = ISAVAILABLE()
            % Check if external dependencies are available for the solver
            bool = true;
        end
        
        function bool = isJavaAvailable()
            % BOOL = ISJAVAAVAILABLE()
            % Check if Java dependencies are available for the solver
            bool = true;
            if ispc % windows
                [~,ret] = dos('java -version');
                if strfind(ret,'not recognized') %#ok<STRIFCND>
                    bool = false;
                end
            else %linux
                [~,ret] = unix('java -version');
                if strfind(ret,'command not found') %#ok<STRIFCND>
                    bool = false;
                end
            end
        end
        
        function fun = accurateStiffOdeSolver()
            % FUN = ACCURATESTIFFODESOLVER()
            % Return default high-accuracy stiff solver
            if isoctave
                %fun = @ode15s;
                fun = @lsode;
            else
                fun = @ode15s;
            end
        end
        
        function fun = accurateOdeSolver()
            % FUN = ACCURATEODESOLVER()
            % Return default high-accuracy non-stiff solver
            if isoctave
                %fun = @ode15s;
                fun = @lsode;
            else
                fun = @ode45;
            end
        end
        
        function fun = fastStiffOdeSolver()
            % FUN = FASTSTIFFODESOLVER()
            % Return default low-accuracy stiff solver
            if isoctave
                %fun = @ode15s;
                fun = @lsode;
            else
                fun = @ode23s;
            end
        end
        
        function fun = fastOdeSolver()
            % FUN = FASTODESOLVER()
            % Return default low-accuracy non-stiff solver
            if isoctave
                %fun = @ode15s;
                fun = @lsode;
            else
                fun = @ode23s;
            end
        end
        
        function [optList, allOpt] = listValidOptions()
            % OPTLIST = LISTVALIDOPTIONS()
            % List valid fields for options data structure
            optList = {'cache','cutoff','force','init_sol','iter_max','iter_tol','tol', ...
                'keep','method','odesolvers','samples','seed','stiff', 'timespan','verbose'};
            allOpt = {'cache','cutoff','force','init_sol','iter_max','iter_tol','tol', ...
                'keep','method','odesolvers','samples','seed','stiff', 'timespan','verbose', ...
                'default','exact','auto','ctmc','ctmc.gpu','gpu','mva','mva.exact','amva','mva.amva',...
                'ssa','ssa.serial.hash','ssa.para.hash','ssa.parallel.hash','ssa.serial','ssa.para',...
                'ssa.parallel','serial.hash','serial','para','parallel','para.hash','parallel.hash',...
                'jmt','jsim','jmva','jmva.mva','jmva.recal','jmva.mom','jmva.comom','jmva.chow','jmva.bs','jmva.aql','jmva.lin','jmva.dmlin','jmva.ls',...
                'jmt.jsim','jmt.jmva','jmt.jmva.mva','jmt.jmva.amva','jmt.jmva.recal','jmt.jmva.comom','jmt.jmva.chow','jmt.jmva.bs','jmt.jmva.aql','jmt.jmva.lin','jmt.jmva.dmlin','jmt.jmva.ls',...
                'ca','comom','mom','propfair','recal','kt', 'nc.comom','nc.gm','nc.mom','nc.propfair','nc.recal','nc.kt', 'nc.ca', ...
                'fluid','nc','nc.exact','nc.imci','nc.ls','nc.le','nc.panacea','panacea','nc.mmint2','mmint2','mam','dec.source','dec.mmap',...
                'mmk','gigk', 'gigk.kingman_approx', ...
                'mm1','mg1','gm1','gig1','gim1','gig1.kingman','gig1.gelenbe','gig1.heyman','gig1.kimura','gig1.allen','gig1.kobayashi','gig1.klb','gig1.marchal','gig1.myskja','gig1.myskja.b',...
                'aba.upper','aba.lower','gb.upper','gb.lower','bjb.upper','bjb.lower','pb.upper','pb.lower'};
        end
        
        function bool = isValidOption(optName)
            % BOOL = ISVALIDOPTION(OPTNAME)
            % Check if the given option exists for the solver
            [~,allOpts] = Solver.listValidOptions();
            bool = any(cell2mat(findstring(optName, allOpts))==1);
        end
        
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            % Return default options
            options = lineDefaults('Solver');
        end
        
        function options = parseOptions(varargin, defaultOptions)
            % OPTIONS = PARSEOPTIONS(VARARGIN, DEFAULTOPTIONS)
            % Parse option parameters into options data structure
            if isempty(varargin)
                options = defaultOptions;
            elseif isstruct(varargin{1})
                options = varargin{1};
            elseif ischar(varargin{1})
                options = defaultOptions;
                [optList, allOpt] = Solver.listValidOptions();
                allMethodsList = setdiff(allOpt, optList);
                while ~isempty(varargin)
                    if Solver.isValidOption(varargin{1})                        
                        switch varargin{1}
                            case allMethodsList
                                options.method = varargin{1};
                                varargin(1) = [];
                            otherwise
                                options.(varargin{1}) = varargin{2};
                                varargin(1) = [];
                                varargin(1) = [];                                
                        end
                    else
                        line_warning(mfilename,'Option "%s" does not exist. Ignoring.',varargin{1});
                        varargin(1) = [];
                    end
                end
            else
                line_error(mfilename,'Invalid parameter.');
            end
        end
        
        function solver = load(chosenmethod, model, varargin)
            % SOLVER = LOAD(CHOSENMETHOD, MODEL, VARARGIN)
            % Returns a solver configured to run the chosen method
            options = Solver.parseOptions(varargin, Solver.defaultOptions);
            options.method = chosenmethod;
            switch options.method
                case {'default','auto'}
                    if strcmp(options.method,'auto'), options.method='default'; end
                    solver = SolverAuto(model, options);
                case {'ctmc','ctmc.gpu','gpu'}
                    if strcmp(options.method,'ctmc'), options.method='default'; end
                    options.method = erase(options.method,'ctmc.');
                    solver = SolverCTMC(model, options);
                case {'mva','mva.exact','amva','mva.amva' ...
                       'mm1','mmk','mg1','gm1','gig1','gim1','gig1.kingman', ...
                       'gigk', 'gigk.kingman_approx', ...
                       'gig1.gelenbe','gig1.heyman','gig1.kimura','gig1.allen','gig1.kobayashi','gig1.klb','gig1.marchal','gig1.myskja','gig1.myskja.b', ...
                       'aba.upper', 'aba.lower', 'bjb.upper', 'bjb.lower', 'gb.upper', 'gb.lower','pb.upper','pb.lower'}
                    if strcmp(options.method,'mva'), options.method='default'; end
                    options.method = erase(options.method,'mva.');
                    solver = SolverMVA(model, options);
                case {'ssa','ssa.serial.hash','ssa.para.hash','ssa.parallel.hash','ssa.serial','ssa.para','ssa.parallel','serial.hash','serial','para','parallel','para.hash','parallel.hash'}
                    if strcmp(options.method,'ssa'), options.method='default'; end
                    options.method = erase(options.method,'ssa.');
                    solver = SolverSSA(model, options);
                case {'jmt','jsim','jmva','jmva.mva','jmva.recal','jmva.comom','jmva.chow','jmva.bs','jmva.aql','jmva.lin','jmva.dmlin','jmva.ls',...
                        'jmt.jsim','jmt.jmva','jmt.jmva.mva','jmt.jmva.amva','jmt.jmva.recal','jmt.jmva.comom','jmt.jmva.chow','jmt.jmva.bs','jmt.jmva.aql','jmt.jmva.lin','jmt.jmva.dmlin','jmt.jmva.ls'}
                    if strcmp(options.method,'jmt'), options.method='default'; end
                    options.method = erase(options.method,'jmt.');
                    solver = SolverJMT(model, options);
                case 'fluid'
                    if strcmp(options.method,'fluid'), options.method='default'; end
                    options.method = erase(options.method,'fluid.');
                    solver = SolverFluid(model, options);
                case {'nc','nc.exact','nc.imci','nc.ls','nc.le','nc.panacea','nc.pana','nc.mmint2','nc.kt','nc.deterministic','nc.sampling','nc.propfair','nc.comom','nc.mom','nc.gm'}
                    if strcmp(options.method,'nc'), options.method='default'; end
                    options.method = erase(options.method,'nc.');
                    solver = SolverNC(model, options);
                case {'mam','mam.dec.source','mam.dec.mmap','mam.dec.poisson'}
                    if strcmp(options.method,'mam'), options.method='default'; end
                    options.method = erase(options.method,'mam.');
                    solver = SolverMAM(model, options);
                otherwise
                    solver = Library(model, options);
            end
        end
    end
end
