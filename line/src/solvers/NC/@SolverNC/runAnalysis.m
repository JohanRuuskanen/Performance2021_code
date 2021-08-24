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
    %line_warning(mfilename,'This model contains features not supported by the solver.');
    ME = MException('Line:FeatureNotSupportedBySolver', 'This model contains features not supported by the solver.');
    throw(ME);
    %                end
    %                runtime = toc(T0);
    %                return
end
Solver.resetRandomGeneratorSeed(options.seed);
[qn] = self.model.getStruct();

if qn.nclosedjobs == 0 && length(qn.nodetype)==3 && all(sort(qn.nodetype)' == sort([NodeType.Source,NodeType.Cache,NodeType.Sink])) % is a non-rentrant cache
    % random initialization
    for ind = 1:qn.nnodes
        if qn.nodetype(ind) == NodeType.Cache
            prob = self.model.nodes{ind}.server.hitClass;
            prob(prob>0) = 0.5;
            self.model.nodes{ind}.server.actualHitProb = prob;
            self.model.nodes{ind}.server.actualMissProb = prob;
        end
    end
    self.model.refreshChains(true);
    % start iteration
    [QN,UN,RN,TN,CN,XN,lG,pij,runtime] = solver_nc_cache_analysis(qn, options);
    self.result.Prob.itemProb = pij;
    for ind = 1:qn.nnodes
        if qn.nodetype(ind) == NodeType.Cache
            %prob = self.model.nodes{ind}.server.hitClass;
            %prob(prob>0) = 0.5;
            for k=1:length(self.model.nodes{ind}.server.hitClass)
                %                for k=1:length(self.model.nodes{ind}.server.hitClass)
                chain_k = find(qn.chains(:,k));
                inchain = find(qn.chains(chain_k,:));
                h = self.model.nodes{ind}.server.hitClass(k);
                m = self.model.nodes{ind}.server.missClass(k);
                if h>0 & m>0
                    self.model.nodes{ind}.server.actualHitProb(k) = XN(h)/nansum(XN(inchain));
                    self.model.nodes{ind}.server.actualMissProb(k) = XN(m)/nansum(XN(inchain));
                end
            end
        end
    end
    self.model.refreshChains;
else % queueing network
    if any(qn.nodetype == NodeType.Cache)
        line_error(mfilename,'Caching analysis not supported yet by NC in general networks.');
    end
    [QN,UN,RN,TN,CN,XN,lG,runtime] = solver_nc_analysis(qn, options);
end
self.setAvgResults(QN,UN,RN,TN,CN,XN,runtime);
self.result.Prob.logNormConstAggr = lG;
end