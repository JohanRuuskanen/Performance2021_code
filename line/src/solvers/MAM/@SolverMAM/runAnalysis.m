function runtime = runAnalysis(self, options, config)
% RUNTIME = RUN()
% Run the solver

T0=tic;
if nargin<2
    options = self.getOptions;
end
if nargin<3
    config = [];
end


if self.enableChecks && ~self.supports(self.model)
    %                if options.verbose
    line_warning(mfilename,'This model contains features not supported by the solver.');
    ME = MException('Line:FeatureNotSupportedBySolver', 'This model contains features not supported by the solver.');
    throw(ME);
    %                end
    %                runtime = toc(T0);
    %                return
end

Solver.resetRandomGeneratorSeed(options.seed);

[qn] = self.model.getStruct();

[Q,U,R,T,C,X] = solver_mam_analysis(qn, options);

runtime=toc(T0);
self.setAvgResults(Q,U,R,T,C,X,runtime);
end