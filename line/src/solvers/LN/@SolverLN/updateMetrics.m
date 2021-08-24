function updateMetrics(self, it)
lqn = self.lqn;

% first obtain svct of activities at hostlayers
self.svct = zeros(self.lqn.nidx,1);
for r=1:size(self.svctmap,1)
    idx = self.svctmap(r,1);
    aidx = self.svctmap(r,2);
    nodeidx = self.svctmap(r,3);
    classidx = self.svctmap(r,4);
    self.svct(aidx) = self.results{end,self.idxhash(idx)}.RN(nodeidx,classidx);
    self.tput(aidx) = self.results{end,self.idxhash(idx)}.TN(nodeidx,classidx);    
    self.svctproc{aidx} = APH.fitMeanAndSCV(self.svct(aidx),3/5);
    %self.svctproc{aidx} = Exp.fitMean(self.svct(aidx));
    %self.tputproc{aidx} = Exp.fitRate(self.tput(aidx));
end

% estimate call response times at hostlayers
self.callrespt = zeros(self.lqn.ncalls,1);
for r=1:size(self.callresptmap,1)
    idx = self.callresptmap(r,1);
    cidx = self.callresptmap(r,2);
    nodeidx = self.callresptmap(r,3);
    classidx = self.callresptmap(r,4);
    self.callrespt(cidx) = self.results{end, self.idxhash(idx)}.RN(nodeidx,classidx);
end

% then resolve the entry svct summing up these contributions
entry_svct = (eye(self.lqn.nidx+self.lqn.ncalls)-self.svctmatrix)\[self.svct;self.callrespt];
entry_svct(1:self.lqn.eshift) = 0;
self.svct(lqn.eshift+1:lqn.eshift+lqn.nentries) = entry_svct(lqn.eshift+1:lqn.eshift+lqn.nentries);
entry_svct((self.lqn.ashift+1):end) = 0;
for r=1:size(self.callresptmap,1)
    cidx = self.callresptmap(r,2);
    eidx = self.lqn.callpair(cidx,2);    
    self.svctproc{eidx} = Exp.fitMean(self.svct(eidx));
end

% determine call response times processes
for r=1:size(self.callresptmap,1)
    %idx = self.callresptmap(r,1);
    cidx = self.callresptmap(r,2);
    %nodeidx = self.callresptmap(r,3);
    %classidx = self.callresptmap(r,4);
    eidx = self.lqn.callpair(cidx,2);
    self.svct(eidx) = entry_svct(eidx);
    self.svctproc{eidx} = Exp.fitMean(self.svct(eidx));
    %ncalls = self.lqn.callproc(cidx).getMean;
    if it==1
        % note that respt is per visit, so number of calls is 1
        self.callrespt(cidx) = self.svct(eidx);
        self.callresptproc{cidx} = self.svctproc{eidx};
    else
        % note that respt is per visit, so number of calls is 1
        self.callresptproc{cidx} = Exp.fitMean(self.callrespt(cidx));
    end
end
end