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

QN = []; UN = [];
RN = []; TN = [];
CN = []; XN = [];
lG = NaN;

if self.enableChecks && ~self.supports(self.model)
   %line_warning(mfilename,'This model contains features not supported by the solver.'); 
ME = MException('Line:FeatureNotSupportedBySolver', 'This model contains features not supported by the solver.'); 
throw(ME);
end
Solver.resetRandomGeneratorSeed(options.seed);

qn = self.model.getStruct();

if (strcmp(options.method,'exact')||strcmp(options.method,'mva')) && ~self.model.hasProductFormSolution
    line_error(mfilename,'The exact method requires the model to have a product-form solution. This model does not have one. You can use Network.hasProductFormSolution() to check before running the solver.');
end

method = options.method;

if qn.nclasses==1 && qn.nclosedjobs == 0 && length(qn.nodetype)==3 && all(sort(qn.nodetype)' == sort([NodeType.Source,NodeType.Queue,NodeType.Sink])) % is a queueing system
    [QN,UN,RN,TN,CN,XN,lG,runtime] = solver_mva_qsys_analysis(qn, options);
elseif qn.nclosedjobs == 0 && length(qn.nodetype)==3 && all(sort(qn.nodetype)' == sort([NodeType.Source,NodeType.Cache,NodeType.Sink])) % is a non-rentrant cache
    % random initialization
    for ind = 1:qn.nnodes
        if qn.nodetype(ind) == NodeType.Cache
            prob = self.model.nodes{ind}.server.hitClass;
            prob(prob>0) = 0.5;
            self.model.nodes{ind}.server.actualHitProb = prob;
            self.model.nodes{ind}.server.actualMissProb = prob;
        end
    end
    self.model.refreshChains;
    % start iteration
    [QN,UN,RN,TN,CN,XN,lG,runtime] = solver_mva_cache_analysis(qn, options);
    
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
            self.model.refreshChains(true);
        end
    end
else % queueing network
    if any(qn.nodetype == NodeType.Cache)
        line_error(mfilename,'Caching analysis not supported yet by MVA in general networks.');
    end
    switch method
        case {'aba.upper', 'aba.lower', 'bjb.upper', 'bjb.lower', 'pb.upper', 'pb.lower', 'gb.upper', 'gb.lower'}
            [QN,UN,RN,TN,CN,XN,lG,runtime] = solver_mva_bound_analysis(qn, options);
        otherwise
            [QN,UN,RN,TN,CN,XN,lG,runtime] = solver_mva_analysis(qn, options);
    end
end
self.setAvgResults(QN,UN,RN,TN,CN,XN,runtime,method);
self.result.Prob.logNormConstAggr = lG;
end
