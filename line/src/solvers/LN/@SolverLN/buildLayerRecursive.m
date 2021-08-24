function  buildLayerRecursive(self, idx, callers, ishostlayer)
lqn = self.lqn;
nreplicas = lqn.repl(idx);
callresptproc = self.callresptproc;
model = Network(lqn.hashnames{idx});
model.setChecks(false); % fast mode
model.attribute = struct('hosts',[],'tasks',[],'entries',[],'activities',[],'calls',[],'serverIdx',0);
if ishostlayer | any(any(lqn.issynccaller(callers, lqn.entriesof{idx}))) %#ok<OR2>
    clientDelay = Delay(model, 'Clients');
    model.attribute.clientIdx = 1;
    model.attribute.serverIdx = 2;
    model.attribute.sourceIdx = NaN;
else
    model.attribute.serverIdx = 1;
    model.attribute.clientIdx = NaN;
    model.attribute.sourceIdx = NaN;
end
serverStation = cell(1,nreplicas);
for m=1:nreplicas
    if m == 1
        serverStation{m} = Queue(model,lqn.hashnames{idx},lqn.sched(idx));
    else
        serverStation{m} = Queue(model,[lqn.hashnames{idx},'.',num2str(m)],lqn.sched(idx));
    end
    serverStation{m}.setNumberOfServers(lqn.mult(idx));
    serverStation{m}.attribute.ishost = ishostlayer;
    serverStation{m}.attribute.idx = idx;
end
aidxClass = cell(1,lqn.nentries+lqn.nacts);
cidxClass = cell(1,0);
cidxAuxClass = cell(1,0);

self.svctmap{idx} = zeros(0,4); % [modelidx, actidx, node, class] % classes to record in svct results at all iterations
self.callresptmap{idx} = zeros(0,4); % [modelidx, callidx, node, class] % call classes to record in respt results at all iterations
self.svcupdmap{idx} = zeros(0,4); % [modelidx, actidx, node, class] % classes to update in the next iteration
self.arvupdmap{idx} = zeros(0,4); % [modelidx, actidx, node, class] % classes to update in the next iteration
self.callupdmap{idx} = zeros(0,4); % [modelidx, callidx, node, class] % calls classes to update in the next iteration
self.routeupdmap{idx} = zeros(0,7); % [modelidx, actidxfrom, actidxto, nodefrom, nodeto, classfrom, classto] % routing probabilities to update in the next iteration

if ishostlayer
    model.attribute.hosts(end+1,:) = [NaN, model.attribute.serverIdx ];
else
    model.attribute.tasks(end+1,:) = [NaN, model.attribute.serverIdx ];
end

hasSource = false; % flag wether a soruce is needed
openClasses = [];
% first pass: create the classes
for tidx_caller = callers
    if ishostlayer | any(any(lqn.issynccaller(tidx_caller, lqn.entriesof{idx}))) %#ok<OR2> % if it is only an asynch caller the closed classes are not needed
        % for each entry of the calling task
        % determine job population
        njobs = lqn.mult(tidx_caller)*lqn.repl(tidx_caller);
        if isinf(njobs)
            callers_of_tidx_caller = find(lqn.taskgraph(:,tidx_caller));
            njobs = sum(lqn.mult(callers_of_tidx_caller)); %#ok<FNDSB>
            if isinf(njobs)
                % if also the callers of tidx_caller are inf servers, then use
                % an heuristic
                njobs = sum(lqn.mult(isfinite(lqn.mult)));
            end
        end
        aidxClass{tidx_caller} = ClosedClass(model, lqn.hashnames{tidx_caller}, njobs, clientDelay);
        aidxClass{tidx_caller}.attribute = [LayeredNetworkElement.TASK, tidx_caller];
        model.attribute.tasks(end+1,:) = [aidxClass{tidx_caller}.index, tidx_caller];
        aidxClass{tidx_caller}.completes = false;
        clientDelay.setService(aidxClass{tidx_caller}, self.thinkproc{tidx_caller});
        if lqn.sched(tidx_caller) ~= SchedStrategy.REF % in 'ref' case the think time is constant
            self.svcupdmap{idx}(end+1,:) = [idx, tidx_caller, 1, aidxClass{tidx_caller}.index];
        end
        for eidx = lqn.entriesof{tidx_caller}
            % create a class
            aidxClass{eidx} = ClosedClass(model, lqn.hashnames{eidx}, 0, clientDelay);
            aidxClass{eidx}.completes = false;
            aidxClass{eidx}.attribute = [LayeredNetworkElement.ENTRY, eidx];
            model.attribute.entries(end+1,:) = [aidxClass{eidx}.index, eidx];
            clientDelay.setService(aidxClass{eidx}, Immediate());
        end
    end
    
    % for each activity of the calling task
    for aidx = lqn.actsof{tidx_caller}
        if ishostlayer | any(any(lqn.issynccaller(tidx_caller, lqn.entriesof{idx}))) %#ok<OR2>
            % create a class
            aidxClass{aidx} = ClosedClass(model, lqn.hashnames{aidx}, 0, clientDelay);
            aidxClass{aidx}.completes = false;
            aidxClass{aidx}.attribute = [LayeredNetworkElement.ACTIVITY, aidx];
            model.attribute.activities(end+1,:) = [aidxClass{aidx}.index, aidx];
            hidx = lqn.parent(lqn.parent(aidx)); % index of host processor
            if ~(ishostlayer && (hidx == idx))
                % set the host demand for the activity
                clientDelay.setService(aidxClass{aidx}, self.svctproc{aidx});
            end
            if ~strcmp(lqn.sched(tidx_caller),SchedStrategy.REF) % in 'ref' case the service activity is constant
                % updmap(end+1,:) = [idx, aidx, 1, idxClass{aidx}.index];
            end
        end
        % add a class for each outgoing call from this activity
        for cidx = lqn.callsof{aidx}
            callmean(cidx) = lqn.callproc{cidx}.getMean;
            switch lqn.calltype(cidx)
                case CallType.ID_ASYNC
                    if lqn.parent(lqn.callpair(cidx,2)) == idx % add only if the target is serverStation
                        if ~hasSource % we need to add source and sink to the model
                            hasSource = true;
                            model.attribute.sourceIdx = length(model.nodes)+1;
                            sourceStation = Source(model,'Source');
                            sinkStation = Sink(model,'Sink');
                        end
                        cidxClass{cidx} = OpenClass(model, lqn.callhashnames{cidx}, 0);
                        sourceStation.setArrival(cidxClass{cidx}, Exp(Distrib.Tol));
                        for m=1:nreplicas
                            serverStation{m}.setService(cidxClass{cidx}, Immediate());
                        end
                        openClasses(end+1,:) = [cidxClass{cidx}.index, callmean(cidx), cidx];
                        model.attribute.calls(end+1,:) = [cidxClass{cidx}.index, cidx, lqn.callpair(cidx,1), lqn.callpair(cidx,2)];
                        aidxClass{cidx}.completes = false;
                        cidxClass{cidx}.attribute = [LayeredNetworkElement.CALL, cidx];
                        minRespT = 0;
                        for tidx_act = lqn.actsof{idx}
                            %self.model.activities{aidx-lqn.ashift}.hostDemand.getMean
                            minRespT = minRespT + lqn.hostdem{tidx_act}.getMean; % upper bound, uses all activities not just the ones reachable by this entry
                        end
                        for m=1:nreplicas
                            serverStation{m}.setService(cidxClass{cidx}, Exp.fitMean(minRespT));
                        end
                    end
                case CallType.ID_SYNC
                    cidxClass{cidx} = ClosedClass(model, lqn.callhashnames{cidx}, 0, clientDelay);
                    model.attribute.calls(end+1,:) = [cidxClass{cidx}.index, cidx, lqn.callpair(cidx,1), lqn.callpair(cidx,2)];
                    aidxClass{cidx}.completes = false;
                    cidxClass{cidx}.attribute = [LayeredNetworkElement.CALL, cidx];
                    minRespT = 0;                    
                    for tidx_act = lqn.actsof{idx}
                        %self.model.activities{aidx-lqn.ashift}.hostDemand.getMean
                        minRespT = minRespT + lqn.hostdem{tidx_act}.getMean; % upper bound, uses all activities not just the ones reachable by this entry
                    end 
                    for m=1:nreplicas
                        serverStation{m}.setService(cidxClass{cidx}, Exp.fitMean(minRespT));
                    end
            end
            
            if callmean(cidx) ~= nreplicas
                switch lqn.calltype(cidx)
                    case CallType.ID_SYNC
                        cidxAuxClass{cidx} = ClosedClass(model, [lqn.callhashnames{cidx},'.Aux'], 0, clientDelay);
                        cidxAuxClass{cidx}.completes = false;
                        cidxAuxClass{cidx}.attribute = [LayeredNetworkElement.CALL, cidx];
                        clientDelay.setService(cidxAuxClass{cidx}, Immediate());
                        for m=1:nreplicas
                            serverStation{m}.setService(cidxAuxClass{cidx}, Disabled());
                        end
                end
            end
        end
    end
end

P = model.initRoutingMatrix;
if hasSource
    for o = 1:size(openClasses,1)
        oidx = openClasses(o,1);
        p = 1 / openClasses(o,2); % divide by mean number of calls, they go to a server at random
        for m=1:nreplicas
            P{model.classes{oidx}, model.classes{oidx}}(sourceStation,serverStation{m}) = 1/nreplicas;
            for n=1:nreplicas
                P{model.classes{oidx}, model.classes{oidx}}(serverStation{m},serverStation{n}) = (1-p)/nreplicas;
            end
            P{model.classes{oidx}, model.classes{oidx}}(serverStation{m},sinkStation) = p;
        end
        cidx = openClasses(o,3);
        self.arvupdmap{idx}(end+1,:) = [idx, cidx, model.getNodeIndex(sourceStation), oidx]; % % 3 = source
        for m=1:nreplicas
            self.callupdmap{idx}(end+1,:) = [idx, cidx, model.getNodeIndex(serverStation{m}), oidx];
            self.callresptmap{idx}(end+1,:) = [idx, cidx, model.getNodeIndex(serverStation{m}), oidx];
        end
    end
end

atClient = 1; % start at client
% second pass: setup the routing out of entries
for tidx_caller = callers
    if lqn.issynccaller(tidx_caller, idx) | ishostlayer % if it is only an asynch caller the closed classes are not needed
        % for each entry of the calling task
        ncaller_entries = length(lqn.entriesof{tidx_caller});
        for eidx = lqn.entriesof{tidx_caller}
            % initialize the probability to select an entry to be identical
            P{aidxClass{tidx_caller}, aidxClass{eidx}}(clientDelay, clientDelay) = 1 / ncaller_entries;
            if ncaller_entries > 1
                % at successive iterations make sure to replace this with throughput ratio
                self.routeupdmap{idx}(end+1,:) = [idx, tidx_caller, eidx, 1, 1, aidxClass{tidx_caller}.index, aidxClass{eidx}.index];
            end
            P = recurActGraph(P, tidx_caller, eidx, aidxClass{eidx}, atClient);
        end
    end
end
model.link(P);
self.ensemble{idx} = model;

    function [P, curClass, jobPos] = recurActGraph(P, tidx_caller, aidx, curClass, jobPos)
        nextaidxs = find(lqn.graph(aidx,:)); % these include the called entries
        %curClass0 = curClass;
        for nextaidx = nextaidxs
            if ~isempty(nextaidx)
                %curClass = curClass0; % for each activity, restart from the input curClass
                if ~(lqn.parent(aidx) == lqn.parent(nextaidx))
                    % if the successor activity is an entry of another task, this is a call
                    cidx = matchrow(lqn.callpair,[aidx,nextaidx]); % find the call index
                    switch lqn.calltype(cidx)
                        case CallType.ID_ASYNC
                            % nop
                        case CallType.ID_SYNC
                            if jobPos == 1
                                if lqn.parent(lqn.callpair(cidx,2)) == idx
                                    % if it is a call to an entry of the server
                                    %P{curClass, curClass}(clientDelay,clientDelay) = full(lqn.graph(lqn.callpair(cidx,1), lqn.callpair(cidx,2))) * (1-pcall);
                                    if callmean(cidx) < nreplicas                                        
                                        P{curClass, cidxAuxClass{cidx}}(clientDelay,clientDelay) = 1 - callmean(cidx);
                                        for m=1:nreplicas
                                            P{curClass, cidxClass{cidx}}(clientDelay,serverStation{m}) = callmean(cidx) / nreplicas;
                                            P{cidxClass{cidx}, cidxClass{cidx}}(serverStation{m},clientDelay) = 1; % not needed, just to avoid leaving the Aux class disconnected
                                        end
                                        P{cidxAuxClass{cidx}, cidxClass{cidx}}(clientDelay,clientDelay) = 1; % not needed, just to avoid leaving the Aux class disconnected
                                    elseif callmean(cidx) == nreplicas
                                        for m=1:nreplicas
                                            P{curClass, cidxClass{cidx}}(clientDelay,serverStation{m}) = 1 / nreplicas;
                                            P{cidxClass{cidx}, cidxClass{cidx}}(serverStation{m},clientDelay) = 1;
                                        end
                                    else % callmean(cidx) > nreplicas
                                        for m=1:nreplicas
                                            P{curClass, cidxClass{cidx}}(clientDelay,serverStation{m}) = 1 / nreplicas;
                                            P{cidxClass{cidx}, cidxAuxClass{cidx}}(serverStation{m},clientDelay) = 1;
                                            P{cidxAuxClass{cidx}, cidxClass{cidx}}(clientDelay,serverStation{m}) = 1 - 1 / (callmean(cidx) / nreplicas);
                                        end
                                        P{cidxAuxClass{cidx}, cidxClass{cidx}}(clientDelay,clientDelay) = 1 / (callmean(cidx));
                                    end
                                    jobPos = 1;
                                    clientDelay.setService(cidxClass{cidx}, Immediate());
                                    for m=1:nreplicas
                                        serverStation{m}.setService(cidxClass{cidx}, callresptproc{cidx});
                                        self.callupdmap{idx}(end+1,:) = [idx, cidx, model.getNodeIndex(serverStation{m}), cidxClass{cidx}.index];
                                        self.callresptmap{idx}(end+1,:) = [idx, cidx, model.getNodeIndex(serverStation{m}), cidxClass{cidx}.index];
                                    end
                                    curClass = cidxClass{cidx};
                                else
                                    % if it is not a call to an entry of the server
                                    if callmean(cidx) < nreplicas
                                        P{curClass, cidxClass{cidx}}(clientDelay,clientDelay) = callmean(cidx)/nreplicas; % the mean number of calls is now embedded in the demand
                                        P{cidxClass{cidx}, cidxAuxClass{cidx}}(clientDelay,clientDelay) = 1; % the mean number of calls is now embedded in the demand
                                        P{curClass, cidxAuxClass{cidx}}(clientDelay,clientDelay) = 1 - callmean(cidx)/nreplicas; % the mean number of calls is now embedded in the demand
                                        curClass = cidxAuxClass{cidx};
                                    elseif callmean(cidx) == nreplicas
                                        P{curClass, cidxClass{cidx}}(clientDelay,clientDelay) = 1;
                                        curClass = cidxClass{cidx};
                                    else % callmean(cidx) > 1
                                        P{curClass, cidxClass{cidx}}(clientDelay,clientDelay) = 1; % the mean number of calls is now embedded in the demand
                                        P{cidxClass{cidx}, cidxClass{cidx}}(clientDelay,clientDelay) = 1 - 1 / (callmean(cidx)/nreplicas); % the mean number of calls is now embedded in the demand
                                        P{cidxClass{cidx}, cidxAuxClass{cidx}}(clientDelay,clientDelay) = 1 / (callmean(cidx)/nreplicas); % the mean number of calls is now embedded in the demand
                                        curClass = cidxAuxClass{cidx};
                                    end
                                    jobPos = 1;
                                    clientDelay.setService(cidxClass{cidx}, callresptproc{cidx});
                                    self.callupdmap{idx}(end+1,:) = [idx, cidx, 1, cidxClass{cidx}.index];
                                    %self.callresptmap{idx}(end+1,:) = [idx, cidx, 1, cidxClass{cidx}.index];
                                end
                            else % job at server
                                if lqn.parent(lqn.callpair(cidx,2)) == idx
                                    % if it is a call to an entry of the server
                                    if callmean(cidx) < nreplicas
                                        for m=1:nreplicas
                                            P{curClass, cidxClass{cidx}}(serverStation{m},clientDelay) = 1 - callmean(cidx);
                                            P{curClass, cidxClass{cidx}}(serverStation{m},serverStation{m}) = callmean(cidx);
                                            serverStation{m}.setService(cidxClass{cidx}, callresptproc{cidx});
                                        end
                                        jobPos = 1;
                                        curClass = cidxAuxClass{cidx};
                                    elseif callmean(cidx) == nreplicas
                                        for m=1:nreplicas
                                            P{curClass, cidxClass{cidx}}(serverStation{m},serverStation{m}) = 1;
                                        end
                                        jobPos = 2;
                                        curClass = cidxClass{cidx};
                                    else % callmean(cidx) > nreplicas
                                        for m=1:nreplicas
                                            P{curClass, cidxClass{cidx}}(serverStation{m},serverStation{m}) = 1;
                                            P{cidxClass{cidx}, cidxClass{cidx}}(serverStation{m},serverStation{m}) = 1 - 1 / (callmean(cidx));
                                            P{cidxClass{cidx}, cidxAuxClass{cidx}}(serverStation{m},clientDelay) = 1 / (callmean(cidx));
                                        end
                                        jobPos = 1;
                                        curClass = cidxAuxClass{cidx};
                                    end
                                    for m=1:nreplicas
                                        serverStation{m}.setService(cidxClass{cidx}, callresptproc{cidx});
                                        self.callupdmap{idx}(end+1,:) = [idx, cidx, model.getNodeIndex(serverStation{m}), cidxClass{cidx}.index];
                                        self.callresptmap{idx}(end+1,:) = [idx, cidx, model.getNodeIndex(serverStation{m}), cidxClass{cidx}.index];
                                    end
                                else
                                    % if it is not a call to an entry of the server
                                    if callmean(cidx) < nreplicas
                                        for m=1:nreplicas
                                            P{curClass, cidxClass{cidx}}(serverStation{m},clientDelay) = callmean(cidx);
                                            P{curClass, cidxAuxClass{cidx}}(serverStation{m},clientDelay) = 1 - callmean(cidx);
                                        end
                                        P{cidxClass{cidx}, cidxAuxClass{cidx}}(clientDelay,clientDelay) = 1;
                                        curClass = cidxAuxClass{cidx};
                                    elseif callmean(cidx) == nreplicas
                                        for m=1:nreplicas
                                            P{curClass, cidxClass{cidx}}(serverStation{m},clientDelay) = 1;
                                        end
                                        curClass = cidxClass{cidx};
                                    else % callmean(cidx) > nreplicas
                                        for m=1:nreplicas
                                            P{curClass, cidxClass{cidx}}(serverStation{m},clientDelay) = 1;
                                        end
                                        P{cidxClass{cidx}, cidxClass{cidx}}(clientDelay,clientDelay) = 1 - 1 / (callmean(cidx));
                                        P{cidxClass{cidx}, cidxAuxClass{cidx}}(clientDelay,clientDelay) = 1 / (callmean(cidx));
                                        curClass = cidxAuxClass{cidx};
                                    end
                                    jobPos = 1;
                                    clientDelay.setService(cidxClass{cidx}, callresptproc{cidx});
                                    self.callupdmap{idx}(end+1,:) = [idx, cidx, 1, cidxClass{cidx}.index];
                                    %self.callresptmap{idx}(end+1,:) = [idx, cidx, 1, cidxClass{cidx}.index];
                                end
                            end
                    end
                end
                % here we have processed all calls, let us do the
                % activities now
                %% if the successor activity is not a call
                if (lqn.parent(aidx) == lqn.parent(nextaidx))
                    if jobPos == 1 % at node 1
                        if ishostlayer
                            for m=1:nreplicas
                                P{curClass, aidxClass{nextaidx}}(clientDelay,serverStation{m}) = full(lqn.graph(aidx,nextaidx));
                                serverStation{m}.setService(aidxClass{nextaidx}, lqn.hostdem{nextaidx});
                            end
                            jobPos = 2;
                            curClass = aidxClass{nextaidx};
                            self.svctmap{idx}(end+1,:) = [idx, nextaidx, 2, aidxClass{nextaidx}.index];
                        else
                            P{curClass, aidxClass{nextaidx}}(clientDelay,clientDelay) = full(lqn.graph(aidx,nextaidx));
                            jobPos = 1;
                            curClass = aidxClass{nextaidx};
                            clientDelay.setService(aidxClass{nextaidx}, self.svctproc{nextaidx});
                            self.svcupdmap{idx}(end+1,:) = [idx, nextaidx, 1, aidxClass{nextaidx}.index];
                        end
                    else % at node 2
                        if ishostlayer
                            for m=1:nreplicas
                                P{curClass, aidxClass{nextaidx}}(serverStation{m},serverStation{m}) = full(lqn.graph(aidx,nextaidx));
                                serverStation{m}.setService(aidxClass{nextaidx}, lqn.hostdem{nextaidx});
                            end
                            jobPos = 2;
                            curClass = aidxClass{nextaidx};
                            self.svctmap{idx}(end+1,:) = [idx, nextaidx, 2, aidxClass{nextaidx}.index];
                        else
                            for m=1:nreplicas
                                P{curClass, aidxClass{nextaidx}}(serverStation{m},clientDelay) = full(lqn.graph(aidx,nextaidx));
                            end
                            jobPos = 1;
                            curClass = aidxClass{nextaidx};
                            clientDelay.setService(aidxClass{nextaidx}, self.svctproc{nextaidx});
                            self.svcupdmap{idx}(end+1,:) = [idx, nextaidx, 1, aidxClass{nextaidx}.index];
                        end
                    end
                    %% now recursively build the rest of the routing matrix graph
                    [P, curClass, jobPos] = recurActGraph(P, tidx_caller, nextaidx, curClass, jobPos);
                    
                    % At this point curClassRec is the last class in the
                    % recursive branch, which we now close with a reply
                    if jobPos == 1
                        P{curClass, aidxClass{tidx_caller}}(clientDelay,clientDelay) = 1;
                        curClass.completes = true;
                    else
                        for m=1:nreplicas
                            P{curClass, aidxClass{tidx_caller}}(serverStation{m},clientDelay) = 1;
                        end
                        curClass.completes = true;
                    end
                end
            end
        end % nextaidx
    end
end
