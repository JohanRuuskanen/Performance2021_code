function nvars = refreshLocalVars(self)
% NVARS = REFRESHLOCALVARS()

nvars = zeros(self.getNumberOfNodes, 1);
varsparam = cell(self.getNumberOfNodes, 1);
rtnodes = self.qn.rtnodes;
for ind=1:self.getNumberOfNodes
    node = self.getNodeByIndex(ind);
    switch class(node)
        case 'Cache'
            nvars(ind) = sum(node.itemLevelCap);
            varsparam{ind} = struct();
            varsparam{ind}.nitems = 0;
            varsparam{ind}.accost = node.accessProb;
            for r=1:self.getNumberOfClasses
                if ~node.popularity{r}.isDisabled
                    varsparam{ind}.nitems = max(varsparam{ind}.nitems,node.popularity{r}.support(2));
                end
            end
            varsparam{ind}.cap = node.itemLevelCap;
            varsparam{ind}.pref = cell(1,self.getNumberOfClasses);
            for r=1:self.getNumberOfClasses
                if node.popularity{r}.isDisabled
                    varsparam{ind}.pref{r} = NaN;
                else
                    varsparam{ind}.pref{r} = node.popularity{r}.evalPMF(1:varsparam{ind}.nitems);
                end
            end
            varsparam{ind}.rpolicy = node.replacementPolicy;
            varsparam{ind}.hitclass = round(node.server.hitClass);
            varsparam{ind}.missclass = round(node.server.missClass);
    end
    switch self.qn.routing(ind)
        case RoutingStrategy.ID_RRB
            nvars(ind) = nvars(ind) + 1;
            % save indexes of outgoing links
            if isempty(varsparam) % reinstantiate if not a cache
                varsparam{ind} = struct();
            end
            varsparam{ind}.outlinks = find(sum(reshape(rtnodes(ind,:)>0,self.qn.nnodes,self.qn.nclasses),2)');
    end
end
if ~isempty(self.qn) %&& isprop(self.qn,'nvars')
    self.qn.setLocalVars(nvars, varsparam);
end
end
