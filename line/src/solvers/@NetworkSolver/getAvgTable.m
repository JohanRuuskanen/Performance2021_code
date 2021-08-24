function [AvgTable,QT,UT,RT,TT] = getAvgTable(self,Q,U,R,T,keepDisabled)
% [AVGTABLE,QT,UT,RT,TT] = GETAVGTABLE(SELF,Q,U,R,T,KEEPDISABLED)
% Return table of average station metrics
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

qn = self.getStruct();

if ~exist('keepDisabled','var')
    keepDisabled = false;
end
M = self.model.getNumberOfStations();
K = self.model.getNumberOfClasses();
if nargin == 2
    keepDisabled = Q;
    [Q,U,R,T,~] = self.model.getAvgHandles();
elseif nargin == 1
    [Q,U,R,T,~] = self.model.getAvgHandles();
end
if isfinite(self.getOptions.timespan(2))
    [Qt,Ut,Tt] = self.model.getTranHandles();
    [QNt,UNt,TNt] = self.getTranAvg(Qt,Ut,Tt);
    QN = cellfun(@(c) c.metric(end),QNt);
    UN = cellfun(@(c) c.metric(end),UNt);
    TN = cellfun(@(c) c.metric(end),TNt);
    RN = zeros(size(QN));
else
    [QN,UN,RN,TN] = self.getAvg(Q,U,R,T);
end

if isempty(QN)
    AvgTable = Table();
    QT = Table();
    UT = Table();
    RT = Table();
    TT = Table();
    AT = Table();
elseif ~keepDisabled
    V = cellsum(qn.visits);
    if isempty(V) % SSA
        for i=1:M
            for c=1:qn.nchains
                chain_classes = find(qn.chains(c,:));
                k = chain_classes(1);
                Tchain=sum(TN(qn.refstat(k),chain_classes));
                for k=chain_classes
                    V(i,k)=TN(i,k)/Tchain;
                end
            end
        end
    end
    
    Qval = []; Uval = [];
    Rval = []; Tval = [];
    Residval =[];
    JobClass = {};
    Station = {};
    for i=1:M
        for k=1:K
            if any(sum([QN(i,k),UN(i,k),RN(i,k),TN(i,k)])>0)
                c = find(qn.chains(:,k));
                JobClass{end+1,1} = Q{i,k}.class.name;
                Station{end+1,1} = Q{i,k}.station.name;
                Qval(end+1) = QN(i,k);
                Uval(end+1) = UN(i,k);
                Rval(end+1) = RN(i,k);
                if RN(i,k)<Distrib.Zero
                    Residval(end+1) = RN(i,k);
                else
                    Residval(end+1) = RN(i,k)*V(i,k)/sum(V(qn.refstat(k),qn.chains(c,:)));
                end
                Tval(end+1) = TN(i,k);
            end
        end
    end
    Station = categorical(Station);
    JobClass = categorical(JobClass);
    QLen = Qval(:); % we need to save first in a variable named like the column
    QT = Table(Station,JobClass,QLen);
    Util = Uval(:); % we need to save first in a variable named like the column
    UT = Table(Station,JobClass,Util);
    RespT = Rval(:); % we need to save first in a variable named like the column
    ResidT = Residval(:); % we need to save first in a variable named like the column
    RT = Table(Station,JobClass,RespT);
    Tput = Tval(:); % we need to save first in a variable named like the column
    TT = Table(Station,JobClass,Tput);
    %AvgTable = Table(Station,JobClass,QLen,Util,RespT,Tput);
    AvgTable = Table(Station, JobClass, QLen, Util, RespT, ResidT, Tput);
else
    V = cellsum(qn.visits);
    if isempty(V) % SSA
        for i=1:M
            for c=1:qn.nchains
                chain_classes = find(qn.chains(c,:));
                k = chain_classes(1);
                Tchain=sum(TN(qn.refstat(k),chain_classes));
                for k=chain_classes
                    V(i,k)=TN(i,k)/Tchain;
                end
            end
        end
    end
    Qval = zeros(M,K); Uval = zeros(M,K);
    Rval = zeros(M,K); Tval = zeros(M,K);
    Residval = zeros(M,K);
    JobClass = cell(K*M,1);
    Station = cell(K*M,1);
    for i=1:M
        for k=1:K
            c = find(qn.chains(:,k));
            JobClass{(i-1)*K+k} = Q{i,k}.class.name;
            Station{(i-1)*K+k} = Q{i,k}.station.name;
            Qval((i-1)*K+k) = QN(i,k);
            Uval((i-1)*K+k) = UN(i,k);
            Rval((i-1)*K+k) = RN(i,k);
            if RN(i,k)<Distrib.Zero
                Residval((i-1)*K+k) = RN(i,k);
            else
                Residval((i-1)*K+k) = RN(i,k)*V(i,k)/sum(V(qn.refstat(k),qn.chains(c,:)));
            end
            Tval((i-1)*K+k) = TN(i,k);
        end
    end
    Station = categorical(Station);
    JobClass = categorical(JobClass);
    QLen = Qval(:); % we need to save first in a variable named like the column
    QT = Table(Station,JobClass,QLen);
    Util = Uval(:); % we need to save first in a variable named like the column
    UT = Table(Station,JobClass,Util);
    RespT = Rval(:); % we need to save first in a variable named like the column
    RT = Table(Station,JobClass,RespT);
    ResidT = Residval(:); % we need to save first in a variable named like the column
    Tput = Tval(:); % we need to save first in a variable named like the column
    TT = Table(Station,JobClass,Tput);
    %AvgTable = Table(Station,JobClass,QLen,Util,RespT,Tput);
    AvgTable = Table(Station, JobClass, QLen, Util, RespT, ResidT, Tput);
end
end
