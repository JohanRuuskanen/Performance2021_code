classdef SolverAuto < NetworkSolver
    %A solver that selects the solution method based on the model characteristics.
    
    %Copyright (c) 2012-2020, Imperial College London
    %All rights reserved.
    
    properties
        CANDIDATE_MAM = 1;
        CANDIDATE_MVA = 2;
        CANDIDATE_NC = 3;
        CANDIDATE_FLUID = 4;
        CANDIDATE_JMT = 5;
        CANDIDATE_SSA = 6;
        CANDIDATE_CTMC = 7;
    end
    
    properties
        candidates; % feasible solvers
        solvers;
    end
    
    methods
        %Constructor
        function self = SolverAuto(model, varargin)
            %SELF = SOLVERAUTO(MODEL, VARARGIN)
            
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            
            %solvers sorted from fastest to slowest
            self.solvers = {};
            self.solvers{1,self.CANDIDATE_MAM} = SolverMAM(model);
            self.solvers{1,self.CANDIDATE_MVA} = SolverMVA(model);
            self.solvers{1,self.CANDIDATE_NC} = SolverNC(model);
            self.solvers{1,self.CANDIDATE_FLUID} = SolverFluid(model);
            self.solvers{1,self.CANDIDATE_JMT} = SolverJMT(model);
            self.solvers{1,self.CANDIDATE_SSA} = SolverSSA(model);
            self.solvers{1,self.CANDIDATE_CTMC} = SolverCTMC(model);
            %turn off warnings temporarily
            wstatus = warning('query');
            warning off;
            boolSolver = [];
            for s=1:length(self.solvers)
                boolSolver(s) = self.solvers{s}.supports(self.model);
                self.solvers{s}.setOptions(self.options);
            end
            self.candidates = {self.solvers{find(boolSolver)}};
            warning(wstatus);
        end
    end
    
    methods
        function bool = supports(self, model)
            %BOOL = SUPPORTS(MODEL)
            
            if isempty(self.candidates)
                bool = false;
            else
                bool = true;
            end
        end
        
        function runtime = runAnalysis(self, options, config) % generic method to run the solver
            %RUNTIME = RUN()
            %Run the solver % GENERIC METHOD TO RUN THE SOLVER
            
            T0 = tic;
            runtime = toc(T0);
            if nargin<2
                options = self.getOptions;
            end
            if nargin<3
                config = [];
            end
            
        end
        
        function [QN,UN,RN,TN] = getAvg(self,Q,U,R,T)
            %[QN,UN,RN,TN] = GETAVG(SELF,Q,U,R,T)
            if nargin ==1
                [Q,U,R,T] = self.model.getAvgHandles;
            elseif nargin == 2
                handlers = Q;
                Q=handlers{1};
                U=handlers{2};
                R=handlers{3};
                T=handlers{4};
            end
            
            % first try with chosen solver, if the method is not available
            %     or fails keep going with the other candidates
            proposedSolvers = {self.chooseSolver(), self.candidates};
            for s=1:length(proposedSolvers)
                try
                    [QN,UN,RN,TN] = proposedSolvers{s}.getAvg(Q,U,R,T);
                    return
                end
            end
        end
        
        function [QNc,UNc,RNc,TNc] = getAvgChain(self,Q,U,R,T)
            %[QNC,UNC,RNC,TNC] = GETAVGCHAIN(SELF,Q,U,R,T)
            
            proposedSolvers = {self.chooseSolver(), self.candidates};
            for s=1:length(proposedSolvers)
                try
                    [QNc,UNc,RNc,TNc] = proposedSolvers{s}.getAvgChain(Q,U,R,T);
                    return
                end
            end
        end
        
        function [CNc,XNc] = getAvgSys(self,R,T)
            %[CNC,XNC] = GETAVGSYS(SELF,R,T)
            
            proposedSolvers = {self.chooseSolver(), self.candidates};
            for s=1:length(proposedSolvers)
                try
                    [CNc,XNc] = proposedSolvers{s}.getAvgSys(R,T);
                    return
                end
            end
        end
        
        function [QNt,UNt,TNt] = getTranAvg(self,Qt,Ut,Tt)
            % [QNT,UNT,TNT] = GETTRANAVG(SELF,QT,UT,TT)
            
            proposedSolvers = {self.chooseSolver(), self.candidates};
            for s=1:length(proposedSolvers)
                try
                    [QNt,UNt,TNt] = proposedSolvers{s}.getTranAvg(Qt,Ut,Tt);
                    return
                end
            end
        end
        
        % AI-based choose solvers
        function solver = chooseSolverAI(self)
            % SOLVER = CHOOSESOLVERAI()
            %
            % This function takes as input a QN model defined in LINE and returns
            % a Solver object with the predicted method loaded
            model = self.model;
            dataVector = extractFeatures(model);
            %Add derived features
            dataVector = [dataVector dataVector(:, 1:3) ./ sum(dataVector(:, 1:3), 2)]; % Percent FCFS, PS, Delay
            dataVector = [dataVector logical(dataVector(:, 4))]; % Has CS or not
            dataVector = [dataVector dataVector(:, 5) ./ sum(dataVector(:, 1:2), 2)]; % Avg svrs per Queue
            dataVector = [dataVector dataVector(:, 7) ./ dataVector(:, 6)]; % Num jobs per chain
            dataVector = [dataVector dataVector(:, 8:10) ./ sum(dataVector(:, 8:10), 2)]; % Percent distributions
            
            load('classifier.mat', 'classifier', 'methodNames', 'selected');
            if isa(classifier, 'cell')
                chosenMethod = predictEnsemble(classifier, dataVector(selected(1:length(dataVector))));
            else
                chosenMethod = predict(classifier, dataVector(selected(1:length(dataVector))));
            end
            
            solver = Solver.load(methodNames(chosenMethod), model);
        end
        
        % chooseSolver: choses a solver from static properties of the model
        function solver = chooseSolver(self)
            % SOLVER = CHOOSESOLVER()
            
            model = self.model;
            if model.hasProductFormSolution
                if model.hasSingleChain
                    %ncoptions = SolverNC.defaultOptions;
                    solver = self.solvers{self.CANDIDATE_NC};
                else % MultiChain
                    if model.hasHomogeneousScheduling(SchedStrategy.INF)
                        solver = self.solvers{self.CANDIDATE_MVA};
                    elseif model.hasMultiServer
                        if sum(model.getNumberOfJobs) / sum(model.getNumberOfChains) > 30 % likely fluid regime
                            solver = self.solvers{self.CANDIDATE_FLUID};
                        elseif sum(model.getNumberOfJobs) / sum(model.getNumberOfChains) > 10 % mid/heavy load
                            solver = self.solvers{self.CANDIDATE_MVA};
                        elseif sum(model.getNumberOfJobs) < 5 % light load
                            solver = self.solvers{self.CANDIDATE_NC};
                        else
                            solver = self.solvers{self.CANDIDATE_MVA};
                        end
                    else % product-form, no infinite servers
                        solver = self.solvers{self.CANDIDATE_NC};
                    end
                end
            else
                solver = chooseSolverAI(self);
                %solver = self.solvers{self.CANDIDATE_MVA};
            end
        end
    end
end