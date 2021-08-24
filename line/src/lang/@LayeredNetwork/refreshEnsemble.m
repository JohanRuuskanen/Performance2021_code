function self = refreshEnsemble(self)
% SELF = REFRESSHENSEMBLE(ISBUILD)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

ensemble = self.ensemble;
lqnGraph = self.getGraph();
graphLayer = self.layerGraph;
clientTask = self.clientTask;

myself = self; % needed for parfor

%% cache recurrent searches
aux = self.aux;
if ~isempty(fieldnames(aux))
    isRefTask = aux.isRefTask;
    isClientHostInSubmodel = aux.isClientHostInSubmodel;
    entries = aux.entries;
    taskobj = aux.taskobj;
    class_entry = aux.class_entry;
    class_hostdemand = aux.class_hostdemand;
else
    isRefTask = cell(length(graphLayer),0);
    isClientHostInSubmodel = cell(length(graphLayer),0);
    taskobj = cell(length(graphLayer),0);
    entries = cell(length(graphLayer),0);
    class_entry = cell(length(graphLayer),0,0);
    class_hostdemand = cell(length(graphLayer),0,0,0);
    for net=1:length(graphLayer)
        serverName = myself.serverName{net};
        jobclass = ensemble{net}.classes;
        for s = 1:length(clientTask{net}) % for all client tasks
            taskobj{net,s} = myself.getNodeObject(clientTask{net}{s});
            isRefTask{net,s} = strcmpi(taskobj{net,s}.scheduling,SchedStrategy.REF);
            isClientHostInSubmodel{net,s} = strcmpi(myself.getNodeHost(clientTask{net}{s}), serverName);
            entries{net,s} = myself.listEntriesOfTask(clientTask{net}{s});
            taskobj{net,s} = myself.getNodeObject(clientTask{net}{s});
            for e=1:length(entries{net,s}) % for all entries{net,s} in the client task
                classEntryName = entries{net,s}{e};
                class_entry{net,s,e} = cellfun(@(c) strcmpi(c.name,classEntryName),jobclass);
                activities = myself.listActivitiesOfEntry(entries{net,s}{e});
                for a=1:length(activities) % for all activities in this entry
                    className = activities{a};                
                    class_hostdemand{net,s,e,a} = cellfun(@(c) strcmpi(c.name,className),jobclass);
                end
            end
        end
    end
    aux.isRefTask = isRefTask;
    aux.isClientHostInSubmodel = isClientHostInSubmodel;
    aux.entries = entries;
    aux.taskobj = taskobj;
    aux.class_entry = class_entry;
    aux.class_hostdemand = class_hostdemand;
    self.aux = aux;
end

%% update parameterization
% this script is valid only for loose layering
for net=1:length(graphLayer)
    hasChanged = false;
    jobclass = ensemble{net}.classes;
    
    %% define stations
    % create surrogate delay node for the layer
    serverName = myself.serverName{net};
    %serverIndex = myself.getNodeIndex(serverName);
    
    %% define routing matrix
    for s = 1:length(clientTask{net}) % for all client tasks
        entries{net,s} = myself.listEntriesOfTask(clientTask{net}{s});
        taskobj{net,s} = myself.getNodeObject(clientTask{net}{s});
        for e=1:length(entries{net,s}) % for all entries{net,s} in the client task
            %entryobj = myself.getNodeObject(entries{net,s}{e});
            %% setup think time for this entry
            classEntryName = entries{net,s}{e};
            if ~isRefTask{net,s}
                eidx = myself.getNodeIndex(classEntryName);
                utilUpperLayerEntry = myself.param.Nodes.Util(eidx);
                tputUpperLayerEntry = myself.param.Nodes.Tput(eidx);
                if utilUpperLayerEntry > 1 % infinite server
                    if utilUpperLayerEntry < 1 + Distrib.Tol
                        utilUpperLayerEntry = 1;
                    else
                        utilUpperLayerEntry = 1;
                        line_warning(mfilename,'Invalid utilization of the upper layer. Setting value to 1.0.');
                    end
                end
                if tputUpperLayerEntry == 0
                    interArrivalFromUpperLayer = Distrib.InfTime;
                else
                    interArrivalFromUpperLayer = jobclass{class_entry{net,s,e}}.population*abs(1-utilUpperLayerEntry) / tputUpperLayerEntry;
                end
                
                thinkTimeMean = taskobj{net,s}.thinkTimeMean;
                %thinkTimeSCV = taskobj{net,s}.thinkTimeSCV;
                destEntryRate = min(Distrib.InfRate,1/(thinkTimeMean + interArrivalFromUpperLayer));
                if thinkTimeMean <= Distrib.Zero
                    destEntryProcess = Exp(destEntryRate);
                    %    [cx,demandMu,demandPhi] = Coxian.fitMeanAndSCV(1/destEntryRate, 1.0);
                elseif interArrivalFromUpperLayer <= Distrib.Zero
                    destEntryProcess = taskobj{net,s}.thinkTime;
                end
                ensemble{net}.nodes{1}.setService(jobclass{class_entry{net,s,e}}, destEntryProcess);
                hasChanged = hasChanged | true;
            end
            
            activities = myself.listActivitiesOfEntry(entries{net,s}{e});
            for a=1:length(activities) % for all activities in this entry
                %% determine properties of this activity
                actobj = myself.getNodeObject(activities{a});
                
                %% first update the host-demand of this activity
                % TODO: spread host-demand in-between calls
                className = activities{a};                
                
                if ~isClientHostInSubmodel{net,s}
                    % if the processor of this client is in another submodel
                    % spend time in the delay equivalent to response time
                    % of this activity
                    destEntryW = myself.param.Nodes.RespT(myself.getNodeIndex(activities{a}));
                    entryRT = Exp.fitMean(destEntryW);
                    ensemble{net}.nodes{1}.setService(jobclass{class_hostdemand{net,s,e,a}}, entryRT);
                    hasChanged = hasChanged | true;
                end
                
                %% update the synchronous calls
                for d=1:length(actobj.syncCallDests)
                    % first return to the delay in the appropriate class
                    destEntry = lqnGraph.Nodes.Name{findstring(lqnGraph.Nodes.Node,actobj.syncCallDests{d})};
                    edgeWeight = myself.edgeWeight(myself.findEdgeIndex(activities{a},destEntry));
                    if edgeWeight >= 1
                        skipProb = 0;
                    else % call-mean < 1
                        skipProb = 1-edgeWeight;
                    end
                    class_synchcall = cellfun(@(c) strcmpi(c.name,[activities{a},'=>',destEntry]),jobclass);
                    
                    % now check if the dest entry's task is server in this submodel
                    isDestTaskInSubmodel = strcmpi(serverName,myself.getNodeTask(destEntry));
                    
                    if isDestTaskInSubmodel
                        % set service time at server for this entry
                        nodeidx = myself.getNodeIndex(destEntry);
                        entryRTObj = Exp.fitMean((1-skipProb)*myself.param.Nodes.RespT(nodeidx));
                        ensemble{net}.nodes{2}.setService(jobclass{class_synchcall}, entryRTObj);
                    else
                        % set as the service time, the response time of this call
                        edgeidx = myself.findEdgeIndex(activities{a},destEntry);
                        destEntryObj = Exp.fitMean((1-skipProb)*myself.param.Edges.RespT(edgeidx));
                        ensemble{net}.nodes{1}.setService(jobclass{class_synchcall}, destEntryObj);
                    end
                    hasChanged = hasChanged | true;
                end
                
                %% update the asynchronous calls
                %for d=1:length(actobj.asyncCallDests)
                    % first return to the delay in the appropriate class
                    %destEntry = lqnGraph.Nodes.Name{findstring(lqnGraph.Nodes.Node,actobj.asyncCallDests{d})};
                    %className = [activities{a},'->',destEntry];
                    %class_asynchcall = cellfun(@(c) strcmpi(c.name,className),jobclass);
                    %destEntryRate = 1/myself.param.Edges.RespT(myself.findEdgeIndex(activities{a},destEntry));
                    % ... to finish
                %end
            end
        end
    end
    if hasChanged
        ensemble{net}.reset();
    end
end
myself.setEnsemble(ensemble);
end

%qn.rates(1,jobclass{class_synchcall}) = destEntryObj.getRate;
%qn.mu{1,jobclass{class_synchcall}} = destEntryObj.getMu;
%qn.phi{1,jobclass{class_synchcall}} = destEntryObj.getPhi;
%qn.proc{1,jobclass{class_synchcall}} = destEntryObj.getRepresentation;
%qn.phases(1,jobclass{class_synchcall}) = length(%qn.phases(1,jobclass{class_synchcall}));
