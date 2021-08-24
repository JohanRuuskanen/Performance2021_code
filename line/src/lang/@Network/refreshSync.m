function sync = refreshSync(self)
% SYNC = REFRESHSYNC()

qn = self.qn;
local = self.getNumberOfNodes+1;
nclasses = qn.nclasses;
sync = {};
emptystate = cellzeros(qn.nnodes,1,0,0);
if any(qn.isstatedep(:))
    rtmask = self.qn.rtfun(emptystate, emptystate);
else
    rtmask = ceil(self.qn.rt);
end

for ind=1:qn.nnodes
    for r=1:qn.nclasses
        if qn.isstation(ind) && qn.phases(qn.nodeToStation(ind),r)> 1
            % Phase-change action
            sync{end+1,1} = struct('active',cell(1),'passive',cell(1));
            sync{end,1}.active{1} = Event(EventType.PHASE, ind, r);
            sync{end,1}.passive{1} = Event(EventType.LOCAL, local, r, 1.0);
        end
        if qn.isstateful(ind)
            if qn.nodetype(ind) == NodeType.Cache
                if ~isnan(qn.varsparam{ind}.pref{r}) % class can read
                    sync{end+1,1}.active{1} = Event(EventType.READ, ind, r);
                    sync{end,1}.passive{1} = Event(EventType.READ, local, r, 1.0);
                end
            end
            isf = qn.nodeToStateful(ind);
            for jnd=1:qn.nnodes
                if qn.isstateful(jnd)
                    jsf = qn.nodeToStateful(jnd);
                    for s=1:nclasses
                        p = rtmask((isf-1)*nclasses+r,(jsf-1)*nclasses+s);
                        if p > 0
                            sync{end+1,1} = struct('active',cell(1),'passive',cell(1));
                            sync{end,1}.active{1} = Event(EventType.DEP, ind, r);
                            switch qn.routing(ind,s)
                                case {RoutingStrategy.ID_RRB, RoutingStrategy.ID_JSQ}
                                    sync{end,1}.passive{1} = Event(EventType.ARV, jnd, s, @(state_before, state_after) at(self.qn.rtfun(state_before, state_after), (isf-1)*nclasses+r, (jsf-1)*nclasses+s));
                                otherwise
                                    sync{end,1}.passive{1} = Event(EventType.ARV, jnd, s, self.qn.rt((isf-1)*nclasses+r, (jsf-1)*nclasses+s));
                            end
                        end
                    end
                end
            end
        end
    end
end
if ~isempty(self.qn) %&& isprop(self.qn,'nvars')
    self.qn.setSync(sync);
end
end
