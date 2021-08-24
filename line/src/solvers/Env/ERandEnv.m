% classdef ERandEnv < ESolver
% % Copyright (c) 2012-2020, Imperial College London
% % All rights reserved.
%
%
%     properties
%         holdingTimes;
%         transitionMatrix;
%         transient_handlers;
%         solversInit;
%         solversIter;
%     end
%
%     methods
%         function self = ERandEnv(networks,holdingTimes,transitionMatrix,options)
% SELF = ERANDENV(NETWORKS,HOLDINGTIMES,TRANSITIONMATRIX,OPTIONS)

%             self@ESolver(networks,options);
%             self.transitionMatrix = transitionMatrix;
%             self.holdingTimes = holdingTimes;
%         end
%
%         function runtime = run(self)
% RUNTIME = RUN()
% Run the solver

%             tic;
%                 options = self.getOptions;
%             if ~isfield(options,'verbose')
%                 options.verbose = 0;
%             end
%
%             if ~isfield(options,'samples')
%                 options.samples = 1e4;
%             end
%
%             if ~isfield(options,'iter_max')
%                 options.iter_max = 100;
%             end
%             if ~isfield(options,'iter_tol')
%                 options.iter_tol = 1e-4;
%             end
%
%             if ~isfield(options,'force')
%                 options.force = false;
%             end
%
%             for e=1:length(self.models)
%                 if ~options.force && ~self.supports(self.models{e})
%                     if options.verbose
%                         line_printf('The model in the random environment stage %d is not supported by the %s solver.',e,self.getName);
%                     end
%                     runtime = toc(T0);
%                     return
%                 end
%             end
%
%             if ~isfield(options,'seed')
%                 options.seed = randi([1,1e6]);
%             end
%
%             if isfield(options,'seed')
%                 rand('seed',options.seed);
%             end
%
%             if ~isfield(options,'timespan')
%                 options.timespan = [0,Inf];
%             end
%
%             E = length(self.models);
%
%             self.solversIter = {};
%             lambda=zeros(1,E);
%             A=zeros(E); I=eye(E);
%             for e=1:E
%                 self.solversInit{e} = SolverMVA(self.models{e},self.options);
%                 self.solversIter{e} = SolverFluid(self.models{e},self.options);
%                 lambda(e) = 1/self.holdingTimes{e}.getMean;
%                 for h=1:E
%                     A(e,h)=-lambda(e)*(I(e,h)-self.transitionMatrix(e,h));
%                 end
%             end
%             pi = ctmc_solve(A);
%             psrc = zeros(E);
%             for e = 1:E
%                 for h = setdiff(1:E,e)
%                     psrc(h,e) = pi(h) * lambda(h) * self.transitionMatrix(h,e);
%                 end
%                 if pi(e) > 0
%                     psrc(:,e) = psrc(:,e) / sum(psrc(:,e));
%                 end
%             end
%
%             for e = 1:E
%                 [Qt,Ut,Xt] = self.models{e}.getTranHandles();
%                 self.models{e}.initFromMarginal(QN);
%                 self.handlers{e} = {Qt,Ut,Xt};
%             end
%
%             %% initialize
%             for it=1:options.iter_max
%                 QE=cell(1,E);
%                 s0 = {};
%                 for e = 1:E
%                     s0(e,:)=self.models{e}.getState;
%                     [QNt{e},UNt{e},TNt{e}] = self.solversIter{e}.getTranAvg(self.transient_handlers{e});
%                     QE{e} = zeros(size(QNt{e}));
%                     for i=1:size(QNt{e},1)
%                         for r=1:size(QNt{e},2)
%                             w{e} = [0; self.holdingTimes{e}.evalCDFInterval(QNt{e}{i,r}.Time(1:end-1), QNt{e}{i,r}.Time(2:end))];
%                             QE{e}(i,r) = QNt{e}{i,r}.Data'*w{e}/sum(w{e});
%                         end
%                     end
%                 end
%                 if it==1
%                     plot(w{e})
%                 end
%                 Q0 = cell(1,E);
%                 for e = 1:E
%                     Q0{e} = zeros(size(QE{e}));
%                     for h=1:E
%                         Q0{e} = Q0{e} + psrc(h,e) * QE{h};
%                     end
%                     self.solversIter{e}.resetResults();
%                     self.models{e}.initFromMarginal(Q0{e});
%                 end
%                 s0 = cell2mat(s0);
%                 if it > 1 && norm(s0-s0_1) < options.iter_tol
%                     break
%                 end
%                 s0_1 = s0;
%             end
%             % compute the other performance indexes
%             for e=1:E
%                 for i=1:size(QNt{e},1)
%                     for r=1:size(QNt{e},2)
%                         w{e} = [0; self.holdingTimes{e}.evalCDFInterval(UNt{e}{i,r}.Time(1:end-1), UNt{e}{i,r}.Time(2:end))];
%                         UE{e}(i,r) = UNt{e}{i,r}.Data'*w{e}/sum(w{e});
%                         w{e} = [0; self.holdingTimes{e}.evalCDFInterval(TNt{e}{i,r}.Time(1:end-1), TNt{e}{i,r}.Time(2:end))];
%                         TE{e}(i,r) = TNt{e}{i,r}.Data'*w{e}/sum(w{e});
%                     end
%                 end
%             end
%
%             QN=0*QE{e};
%             UN=0*UE{e};
%             TN=0*TE{e};
%             for e=1:E
%                 QN = QN + pi(e) * QE{e};
%                 UN = UN + pi(e) * UE{e};
%                 TN = TN + pi(e) * TE{e};
%             end
%             runtime = toc;
%             self.result.('solver') = self.getName();
%             self.result.Avg.Q = QN;
%             %    self.result.Avg.R = R;
%             %    self.result.Avg.X = X;
%             self.result.Avg.U = UN;
%             self.result.Avg.T = TN;
%             %    self.result.Avg.C = C;
%             self.result.runtime = runtime;
%         end
%
%         function name = getName(self)
% NAME = GETNAME()

%             name = mfilename;
%         end
%
%         function [results,runtime] = solve(self, options)
% [RESULTS,RUNTIME] = SOLVE(OPTIONS)

%             %            try
%                 options = self.getOptions;
%             runtime = self.runAnalysis();
%             results = self.getResults();
%             %            catch me
%             %                line_error(getReport(me));
%             %            end
%         end
%
%         function [QNclass, UNclass, TNclass] = getAvg(self)
% [QNCLASS, UNCLASS, TNCLASS] = GETAVG()

%             if isempty(self.result) || (isfield(self.options,'force') && self.options.force)
%                 self.solve(self.options);
%                 if isempty(self.result)
%                     QNclass=[];
%                     UNclass=[];
%                     TNclass=[];
%                     return
%                 end
%             end
%             QNclass = self.result.Avg.Q;
%             UNclass = self.result.Avg.U;
%             TNclass = self.result.Avg.T;
%         end
%     end
%
%     methods (Static)
%         function [bool, featSupported] = supports(model)
% [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)

%             featUsed = model.getUsedLangFeatures();
%
%             featSupported = SolverFeatureSet;
%
%             % Nodes
%             featSupported.setTrue('ClassSwitch');
%             featSupported.setTrue('DelayStation');
%             featSupported.setTrue('Queue');
%             featSupported.setTrue('Sink');
%             featSupported.setTrue('Source');
%
%             % Distributions
%             featSupported.setTrue('Cox2');
%             featSupported.setTrue('Erlang');
%             featSupported.setTrue('Exponential');
%             featSupported.setTrue('HyperExp');
%
%             % Sections
%             featSupported.setTrue('StatelessClassSwitcher'); % Section
%             featSupported.setTrue('InfiniteServer'); % Section
%             featSupported.setTrue('SharedServer'); % Section
%             featSupported.setTrue('Buffer'); % Section
%             featSupported.setTrue('Dispatcher'); % Section
%             featSupported.setTrue('Server'); % Section (Non-preemptive)
%             featSupported.setTrue('JobSink'); % Section
%             featSupported.setTrue('RandomSource'); % Section
%             featSupported.setTrue('ServiceTunnel'); % Section
%
%             % Scheduling strategy
%             featSupported.setTrue('SchedStrategy_INF');
%             featSupported.setTrue('SchedStrategy_PS');
%
%             % Customer Classes
%             featSupported.setTrue('ClosedClass');
%             featSupported.setTrue('OpenClass');
%
%             bool = SolverFeatureSet.supports(featSupported, featUsed);
%         end
%     end
% end
%
