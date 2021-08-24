classdef SolverJMT < NetworkSolver
    % A solver that interfaces the Java Modelling Tools (JMT) to LINE.
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    %Private properties
    properties %(GetAccess = 'private', SetAccess='private')
        jmtPath;
        filePath;
        fileName;
        maxSimulatedTime;
        maxSamples;
        maxEvents;
        seed;
    end
    
    %Constants
    properties (Constant)
        xmlnsXsi = 'http://www.w3.org/2001/XMLSchema-instance';
        xsiNoNamespaceSchemaLocation = 'Archive.xsd';
        fileFormat = 'jsimg';
        jsimgPath = '';
    end
    
    % PUBLIC METHODS
    methods
        
        %Constructor
        function self = SolverJMT(model, varargin)
            % SELF = SOLVERJMT(MODEL, VARARGIN)
            
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            if ~Solver.isJavaAvailable
                line_error(mfilename,'SolverJMT requires the java command to be available on the system path.');
            end
            if ~Solver.isAvailable
                line_error(mfilename,'SolverJMT cannot located JMT.jar in the MATLAB path.');
            end
            
            self.maxEvents = -1;
            jarPath = jmtGetPath;
            self.setJMTJarPath(jarPath);
            filePath = tempdir;
            self.filePath = filePath;
            [~,fileName]=fileparts(tempname);
            self.fileName = fileName;
        end
        
        [simDoc, section] = saveArrivalStrategy(self, simDoc, section, currentNode)
        [simDoc, section] = saveBufferCapacity(self, simDoc, section, currentNode)
        [simDoc, section] = saveClassSwitchStrategy(self, simDoc, section, currentNode)
        [simDoc, section] = saveDropStrategy(self, simDoc, section)
        [simDoc, section] = saveForkStrategy(self, simDoc, section, currentNode)
        [simDoc, section] = saveGetStrategy(self, simDoc, section, currentNode)
        [simDoc, section] = saveJoinStrategy(self, simDoc, section, currentNode)
        [simDoc, section] = saveLogTunnel(self, simDoc, section, currentNode)
        [simDoc, section] = saveNumberOfServers(self, simDoc, section, currentNode)
        [simDoc, section] = savePreemptiveStrategy(self, simDoc, section, currentNode)
        [simDoc, section] = savePreemptiveWeights(self, simDoc, section, currentNode)
        [simDoc, section] = savePutStrategy(self, simDoc, section, currentNode)
        [simDoc, section] = saveRoutingStrategy(self, simDoc, section, currentNode)
        [simDoc, section] = saveServerVisits(self, simDoc, section)
        [simDoc, section] = saveServiceStrategy(self, simDoc, section, currentNode)
        [simElem, simDoc] = saveClasses(self, simElem, simDoc)
        [simElem, simDoc] = saveLinks(self, simElem, simDoc)
        [simElem, simDoc] = saveMetrics(self, simElem, simDoc)
        [simElem, simDoc] = saveXMLHeader(self, logPath)
        
        function fileName = getFileName(self)
            % FILENAME = GETFILENAME()
            
            fileName = self.fileName;
        end
        
        %Setter
        function self = setJMTJarPath(self, path)
            % SELF = SETJMTJARPATH(PATH)
            
            self.jmtPath = path;
        end
        
        % Getters
        function out = getJMTJarPath(self)
            % OUT = GETJMTJARPATH()
            
            out = self.jmtPath;
        end
        
        function out = getFilePath(self)
            % OUT = GETFILEPATH()
            
            out = self.filePath;
        end
        
        Tsim = run(self, options)
        
        jwatView(self, options)
        jsimgView(self, options)
        
        [outputFileName] = writeJMVA(self, outputFileName)
        [outputFileName] = writeJSIM(self, outputFileName)
        [result, parsed] = getResults(self)
        [result, parsed] = getResultsJSIM(self)
        [result, parsed] = getResultsJMVA(self)
    end
    
    %Private methods.
    methods (Access = 'private')
        function out = getJSIMTempPath(self)
            % OUT = GETJSIMTEMPPATH()
            
            fname = [self.getFileName(), ['.', 'jsimg']];
            out = [self.filePath,'jsimg',filesep, fname];
        end
        
        function out = getJMVATempPath(self)
            % OUT = GETJMVATEMPPATH()
            
            fname = [self.getFileName(), ['.', 'jmva']];
            out = [self.filePath,'jmva',filesep, fname];
        end
    end
    
    %Private methods.
    methods (Access = 'protected')
        function bool = hasAvgResults(self)
            % BOOL = HASAVGRESULTS()
            
            bool = self.hasResults();
        end
    end
    
    
    methods (Access = 'public')
        getProbNormConstAggr(self); % jmva
        %% StateAggr methods
        Pr = getProbAggr(self, node, state_a);
        [Pi_t, SSnode_a] = getTranProbAggr(self, node);
        probSysStateAggr = getProbSysAggr(self);
        tranNodeStateAggr = sampleAggr(self, node, numsamples);
        tranSysStateAggr = sampleSysAggr(self, numsamples);
        
        %% Cdf methods
        RD = getCdfRespT(self, R);
        RD = getTranCdfRespT(self, R);
        RD = getTranCdfPassT(self, R);
    end
    
    
    methods (Static)
        
        function bool = isAvailable()
            % BOOL = ISAVAILABLE()
            
            bool = true;
            if isempty(which('JMT.jar'))
                bool = false;
            end
        end
        
        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()
            
            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink',...
                'Source',...
                'Router',...
                'ClassSwitch',...
                'DelayStation',...
                'Queue',...
                'Fork',...
                'Join',...
                'Forker',...
                'Joiner',...
                'Logger',...
                'Coxian',...
                'Cox2',...
                'APH',...
                'Erlang',...
                'Exponential',...
                'HyperExp',...
                'Det',...
                'Gamma',...
                'MAP',...
                'MMPP2',...
                'Normal',...
                'Pareto',...
                'Replayer',...
                'Uniform',...
                'StatelessClassSwitcher',...
                'InfiniteServer',...
                'SharedServer',...
                'Buffer',...
                'Dispatcher',...
                'Server',...
                'JobSink',...
                'RandomSource',...
                'ServiceTunnel',...
                'LogTunnel',...
                'SchedStrategy_INF',...
                'SchedStrategy_PS',...
                'SchedStrategy_DPS',...
                'SchedStrategy_FCFS',...
                'SchedStrategy_GPS',...
                'SchedStrategy_RAND',...
                'SchedStrategy_HOL',...
                'SchedStrategy_LCFS',...
                'SchedStrategy_SEPT',...
                'SchedStrategy_LEPT',...
                'SchedStrategy_SJF',...
                'SchedStrategy_LJF',...
                'RoutingStrategy_PROB',...
                'RoutingStrategy_RAND',...
                'RoutingStrategy_RRB',...
                'SchedStrategy_EXT',...
                'ClosedClass',...
                'OpenClass'});
        end
        
        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverJMT.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end
        
        function jsimgOpen(filename)
            % JSIMGOPEN(FILENAME)
            
            [path] = fileparts(filename);
            if isempty(path)
                filename=[pwd,filesep,filename];
            end
            runtime = java.lang.Runtime.getRuntime();
            cmd = ['java -cp "',jmtGetPath,filesep,'JMT.jar" jmt.commandline.Jmt jsimg "',filename,'"'];
            system(cmd);
            %runtime.exec(cmd);
        end
        
        function jsimwOpen(filename)
            % JSIMWOPEN(FILENAME)
            
            runtime = java.lang.Runtime.getRuntime();
            cmd = ['java -cp "',jmtGetPath,filesep,'JMT.jar" jmt.commandline.Jmt jsimw "',which(filename)];
            %system(cmd);
            runtime.exec(cmd);
        end
        
        dataSet = parseLogs(model, isNodeLogged, metric);
        [state, evtype, evclass] = parseTranState(fileArv, fileDep, nodePreload);
        [classResT, jobResT, jobResTArvTS, classResTJobID] = parseTranRespT(fileArv, fileDep);
        
        function checkOptions(options)
            % CHECKOPTIONS(OPTIONS)
            
            solverName = mfilename;
            %             if isfield(options,'timespan')  && isfinite(options.timespan(2))
            %                 line_error(mfilename,'Finite timespan not supported in %s',solverName);
            %             end
        end
        
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = lineDefaults('JMT');
        end
    end
    
end
