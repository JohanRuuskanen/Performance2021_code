function updateLayers(self, it)
lqn = self.lqn;
ensemble = self.ensemble;
idxhash = self.idxhash;
idletproc = self.idletproc;
svctproc = self.svctproc;
svcupdmap = self.svcupdmap;
arvupdmap = self.arvupdmap;
tputproc = self.tputproc;
callupdmap = self.callupdmap;
callresptproc = self.callresptproc;

% reassign svctimes
for r=1:size(svcupdmap,1)
    if mod(it, 0)
        ri = size(svcupdmap,1) - r + 1;
    else
        ri = r;
    end
    idx = svcupdmap(ri,1);
    aidx = svcupdmap(ri,2);
    nodeidx = svcupdmap(ri,3);
    classidx = svcupdmap(ri,4);
    class = ensemble{idxhash(svcupdmap(ri,1))}.classes{classidx};
    switch nodeidx
        case  ensemble{idxhash(idx)}.attribute.clientIdx
            if lqn.type(aidx) == LayeredNetworkElement.TASK
                if lqn.sched(aidx) ~= SchedStrategy.REF
                    if ~isempty(idletproc{aidx}) % this is empty for isolated components, which can be ignored
                        ensemble{idxhash(idx)}.nodes{nodeidx}.setService(class, idletproc{aidx});
                    end
                else
                    ensemble{idxhash(idx)}.nodes{nodeidx}.setService(class, svctproc{aidx});
                end
            else
                ensemble{idxhash(idx)}.nodes{nodeidx}.setService(class, svctproc{aidx});
            end
        case ensemble{idxhash(idx)}.attribute.serverIdx
            ensemble{idxhash(idx)}.nodes{nodeidx}.setService(class, svctproc{aidx});
    end
end

% reassign arvrates
for r=1:size(arvupdmap,1)
    if mod(it, 0)
        ri = size(arvupdmap,1) - r + 1;
    else
        ri = r;
    end
    idx = arvupdmap(ri,1);
    cidx = arvupdmap(ri,2);
    nodeidx = arvupdmap(ri,3);
    classidx = arvupdmap(ri,4);
    class = ensemble{idxhash(arvupdmap(ri,1))}.classes{classidx};
    ensemble{idxhash(idx)}.nodes{nodeidx}.setArrival(class, tputproc{lqn.callpair(cidx,1)});
end

% reassign call svct / respt
for c=1:size(callupdmap,1)
    if mod(it, 0)
        ci = size(callupdmap,1) - c + 1;
    else
        ci = c;
    end
    idx = callupdmap(ci,1);
    cidx = callupdmap(ci,2);
    nodeidx = callupdmap(ci,3);
    class = ensemble{idxhash(callupdmap(ci,1))}.classes{callupdmap(ci,4)};
    switch nodeidx
        case ensemble{idxhash(idx)}.attribute.clientIdx % client
             ensemble{idxhash(idx)}.nodes{nodeidx}.setService(class, callresptproc{cidx});
        case ensemble{idxhash(idx)}.attribute.serverIdx % the call is processed by the server, then replace with the svc time
            eidx = lqn.callpair(cidx,2);
            ensemble{idxhash(idx)}.nodes{nodeidx}.setService(class, svctproc{eidx});
    end
end


%self.ensemble = ensemble;
%self.idxhash = idxhash;
%self.idletproc = idletproc;
%self.svctproc = svctproc;
%self.svcupdmap = svcupdmap;
%self.arvupdmap = arvupdmap;
%self.tputproc = tputproc;
%self.callupdmap = callupdmap;
%self.callresptproc = callresptproc;
end