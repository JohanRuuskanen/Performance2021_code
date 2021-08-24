function [runtime, tranSysState, tranSync] = run(self, options)
% [RUNTIME, TRANSYSSTATE] = RUN()

T0=tic;
if ~exist('options','var')
    options = self.getOptions;
end

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

qn = self.getStruct();

% TODO: add priors on initial state
qn.state = self.model.getState; % not used internally by SSA

[Q,U,R,T,C,X,~, tranSysState, tranSync, qn] = solver_ssa_analysis(qn, options);
for isf=1:qn.nstateful   
    ind = qn.statefulToNode(isf);
    switch class(self.model.nodes{qn.statefulToNode(isf)})
        case 'Cache'
            self.model.nodes{qn.statefulToNode(isf)}.server.actualHitProb = qn.varsparam{ind}.actualhitprob;
            self.model.nodes{qn.statefulToNode(isf)}.server.actualMissProb = qn.varsparam{ind}.actualmissprob;
            self.model.refreshChains(true);
    end
end
runtime = toc(T0);
self.setAvgResults(Q,U,R,T,C,X,runtime);
self.result.space = qn.space;
end