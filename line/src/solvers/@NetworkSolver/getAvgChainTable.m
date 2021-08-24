function [AvgChain,QTc,UTc,RTc,TTc] = getAvgChainTable(self,Q,U,R,T)
% [AVGCHAIN,QTC,UTC,RTC,TTC] = GETAVGCHAINTABLE(SELF,Q,U,R,T)

% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

M = self.model.getNumberOfStations();
K = self.model.getNumberOfChains();
if nargin == 1
    [Q,U,R,T] = self.model.getAvgHandles();
end
[QNc,UNc,RNc,TNc] = self.getAvgChain(Q,U,R,T);

ChainObj = self.model.getChains();
ChainName = cellfun(@(c) c.name,ChainObj,'UniformOutput',false);
ChainClasses = cell(1,length(ChainName));
for c=1:length(ChainName)
    ChainClasses{c} = ChainObj{c}.classnames;
end
if isempty(QNc)
    AvgChain = Table();
    QTc = Table();
    UTc = Table();
    RTc = Table();
    TTc = Table();
else
    Qval = zeros(M,K); Uval = zeros(M,K);
    Rval = zeros(M,K); Tval = zeros(M,K);
    Chain = cell(K*M,1);
    JobClasses = cell(K*M,1);
    Station = cell(K*M,1);
    for i=1:M
        for k=1:K
            Chain{(i-1)*K+k} = ChainName{k};
            JobClasses((i-1)*K+k,1) = {categorical(ChainClasses{k}(:))};
            Station{(i-1)*K+k} = Q{i,k}.station.name;
            Qval((i-1)*K+k) = QNc(i,k);
            Uval((i-1)*K+k) = UNc(i,k);
            Rval((i-1)*K+k) = RNc(i,k);
            Tval((i-1)*K+k) = TNc(i,k);
        end
    end
    Chain = categorical(Chain);
    Station = categorical(Station);
    QLen = Qval(:); % we need to save first in a variable named like the column
    QTc = Table(Station,Chain,JobClasses,QLen);
    Util = Uval(:); % we need to save first in a variable named like the column
    UTc = Table(Station,Chain,JobClasses,Util);
    RespT = Rval(:); % we need to save first in a variable named like the column
    RTc = Table(Station,Chain,JobClasses,RespT);
    Tput = Tval(:); % we need to save first in a variable named like the column
    TTc = Table(Station,Chain,JobClasses,Tput);
    AvgChain = Table(Station,Chain,JobClasses,QLen,Util,RespT,Tput);
end
end
