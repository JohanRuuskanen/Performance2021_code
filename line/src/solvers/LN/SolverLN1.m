classdef SolverLN1 < LayeredNetworkSolver & EnsembleSolver
    % LINE native solver for layered networks.
    %
    % Copyright (c) 2012-2020, Imperial College London
    % All rights reserved.
    
    properties
        initNodeAvgTable = {};
        initCallAvgTable = {};
    end
    
    properties (Hidden)
        successors = {};
    end
    
    methods
        function self = reset(self)
            self.initNodeAvgTable = {};
            self.initCallAvgTable = {};
        end
        
        function self = SolverLN1(model, solverFactory, varargin)
            % SELF = SOLVERLN1(MODEL,SOLVERFACTORY,VARARGIN)
            
            self@LayeredNetworkSolver(model, mfilename);
            self@EnsembleSolver(model, mfilename);
            
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            self.model.initDefault;
            self.ensemble = self.model.getEnsemble();
            for e=1:self.getNumberOfModels
                if nargin == 1
                    self.setSolver(SolverMVA(self.ensemble{e}),e);
                else
                    self.setSolver(solverFactory(self.ensemble{e}),e);
                end
            end
            %            model.refreshLayers;
        end
        
        function initFromRawAvgTables(self, NodeAvgTable, CallAvgTable)
            self.initNodeAvgTable = NodeAvgTable;
            self.initCallAvgTable = CallAvgTable;
        end
        
        function paramFromRawAvgTables(self, NodeAvgTable, CallAvgTable)
            self.model.param.Nodes.RespT = NodeAvgTable.Phase1ServiceTime;
            self.model.param.Nodes.RespT(isnan(self.model.param.Nodes.RespT)) = 0;
            
            self.model.param.Nodes.Tput = NodeAvgTable.Throughput;
            self.model.param.Nodes.Tput(isnan(self.model.param.Nodes.Tput)) = 0;
            for e = find(strcmpi(NodeAvgTable.NodeType,'Entry'))
                self.model.param.Nodes.Util(e) = NodeAvgTable.Phase1Utilization(e);
            end
            
            for h = 1:length(CallAvgTable.SourceNode)
                edgeidx = self.model.findEdgeIndex(CallAvgTable.SourceNode{h}, CallAvgTable.TargetNode{h});
                etargetidx = self.model.getNodeIndex(CallAvgTable.TargetNode{h});
                self.model.param.Edges.RespT(edgeidx) = CallAvgTable.Waiting(h) + self.model.param.Nodes.RespT(etargetidx);
                self.model.param.Nodes.RespT(edgeidx) = abs(self.model.param.Nodes.RespT(edgeidx) - self.model.param.Edges.RespT(edgeidx));
            end
        end
        
        function bool = converged(self, it) % convergence test at iteration it
            % BOOL = CONVERGED(IT) % CONVERGENCE TEST AT ITERATION IT
            
            bool = false;
            if it > 1
                maxIterErr = 0;
                for e = 1:size(self.results,2)
                    try
                        maxIterErr = max([maxIterErr, nanmax(abs(1 - self.results{end,e}.RespT ./ self.results{end-1,e}.RespT))]);
                    catch
                        maxIterErr = Inf;
                        break
                    end
                end
                if self.options.verbose == 2
                    line_printf( sprintf('SolverLN error is: %f',maxIterErr));
                end
                if maxIterErr < self.options.iter_tol
                    if self.options.verbose
                        line_printf( sprintf('\nSolverLN completed in %d iterations.',size(self.results,1)));
                    end
                    bool = true;
                end
            end
        end
        
        function init(self) % operations before starting to iterate
            % INIT() % OPERATIONS BEFORE STARTING TO ITERATE
            %nop
        end
        
        function pre(self, it) % operations before an iteration
            % PRE(IT) % OPERATIONS BEFORE AN ITERATION
            
            %nop
        end
        
        function [result, runtime] = analyze(self, it, e)
            % [RESULT, RUNTIME] = ANALYZE(IT, E)
            T0 = tic;
            %self.ensemble{e}.reset();
            self.solvers{e}.resetResults();
            if it>1  
                self.ensemble{e}.initFromAvgTableQLen(self.results{it-1,e});
            end
            result = self.solvers{e}.getAvgTable(true); % return also zero metrics
            runtime = toc(T0);
        end
        
        function post(self, it) % operations after an iteration
            % POST(IT) % OPERATIONS AFTER AN ITERATION
            for netSortAscending = [false, true] % do elevator up and down
                self.model.updateParam({self.results{it,:}}, netSortAscending);
            end
            if it == 1 && ~isempty(self.initNodeAvgTable)
                self.paramFromRawAvgTables(self.initNodeAvgTable, self.initCallAvgTable);
            end
            self.ensemble = self.model.refreshLayers(); % update Network objects in ensemble
        end
        
        function finish(self) % operations after iterations are completed
            % FINISH() % OPERATIONS AFTER INTERATIONS ARE COMPLETED
            if self.options.verbose
                line_printf('\n');
            end            
        end
        
        function [QN,UN,RN,TN] = getAvg(self,~,~,~,~)
            % [QN,UN,RN,TN] = GETAVG(SELF,~,~,~,~)
            
            self.runAnalysis(); % run iterations
            lqnGraph = self.model.getGraph;
            Avg = self.model.param;
            % At this point the Avg data structure includes only the
            % fundamental perf indexes that uniquely determine a valid LQN.
            % We now derive the other perf indexes
            for edge = 1:height(lqnGraph.Edges)
                if lqnGraph.Edges.Type(edge) == 1 % add contribution of sync-calls
                    syncSource = lqnGraph.Edges.EndNodes{edge,1};
                    aidx = findstring(lqnGraph.Nodes.Name,syncSource);
                    if lqnGraph.Edges.Weight(edge) >= 1
                        Avg.Nodes.RespT(aidx) = Avg.Nodes.RespT(aidx) +  Avg.Edges.RespT(edge) * lqnGraph.Edges.Weight(edge);
                    else
                        Avg.Nodes.RespT(aidx) = Avg.Nodes.RespT(aidx) +  Avg.Edges.RespT(edge);
                    end
                end
            end
            % - qlen is respT * tput
            % - qlen of task is sum of qlen of its entries
            % - tput of task is sum of tput of its entries
            % - util of entry is sum of util of its activities
            Avg.Nodes.QLen = Avg.Nodes.RespT .* Avg.Nodes.Tput;
            Avg.Nodes.Util = lqnGraph.Nodes.D .* Avg.Nodes.Tput;
            procPos = strcmp(lqnGraph.Nodes.Type,'H');
            Avg.Nodes.QLen(procPos) = NaN;
            Avg.Nodes.RespT(procPos) = NaN;
            Avg.Nodes.Tput(procPos) = NaN;
            taskPos = strcmp(lqnGraph.Nodes.Type,'R') | strcmp(lqnGraph.Nodes.Type,'T');
            for tidx = find(taskPos)'
                entriesOfTask = self.model.listEntriesOfTask(tidx);
                for e = 1:length(entriesOfTask)
                    Avg.Nodes.Tput(tidx) = 0;
                end
                for e = 1:length(entriesOfTask)
                    eidx = self.model.getNodeIndex(entriesOfTask{e});
                    Avg.Nodes.QLen(tidx) = Avg.Nodes.QLen(tidx) + Avg.Nodes.QLen(eidx);
                    Avg.Nodes.Tput(tidx) = Avg.Nodes.Tput(tidx) + Avg.Nodes.Tput(eidx);
                    actOfEntry = self.model.listActivitiesOfEntry(entriesOfTask{e});
                    Avg.Nodes.Util(eidx) = 0;
                    for a = 1:length(actOfEntry)
                        aidx = self.model.getNodeIndex(actOfEntry{a});
                        Avg.Nodes.Util(eidx) = Avg.Nodes.Util(eidx) + Avg.Nodes.Util(aidx);
                    end
                    Avg.Nodes.Util(tidx) = Avg.Nodes.Util(tidx) + Avg.Nodes.Util(eidx);
                end
                pidx = self.model.getNodeIndex(lqnGraph.Nodes.Proc{tidx});
                Avg.Nodes.Util(pidx) = Avg.Nodes.Util(pidx) + Avg.Nodes.Util(tidx);
            end
            Avg.Nodes.RespT(taskPos) = NaN;
            QN =  Avg.Nodes.QLen;
            UN =  Avg.Nodes.Util;
            TN =  Avg.Nodes.Tput;
            RN =  Avg.Nodes.RespT;
        end
    end
    
    methods (Static)
        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink','Source','Queue',...
                'Coxian','Erlang','Exponential','HyperExp',...
                'Buffer','Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_PS','SchedStrategy_FCFS','ClosedClass'});
            bool = true;
            for e = 1:model.getNumberOfLayers()
                bool = bool && SolverFeatureSet.supports(featSupported, featUsed{e});
            end
        end
    end
    
    methods (Static)
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = lineDefaults('LN');
        end
    end
end
