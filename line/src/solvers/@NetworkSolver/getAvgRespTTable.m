function [AvgTable,T] = getAvgRespTTable(self,R,keepDisabled)
% [AVGTABLE,T] = GETAVGRESPTTABLE(SELF,R,KEEPDISABLED)

% Return table of average station metrics
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
if ~exist('keepDisabled','var')
    keepDisabled = false;
end

M = self.model.getNumberOfStations();
K = self.model.getNumberOfClasses();
if nargin == 1
    R = self.model.getAvgRespTHandles();
end
RN = self.getAvgRespT();
if isempty(RN)
    AvgTable = Table();
    RT = Table();
elseif ~keepDisabled
    Rval = [];
    JobClass = {};
    Station = {};
    for i=1:M
        for k=1:K
            if any(sum([RN(i,k)])>0)
                JobClass{end+1,1} = R{i,k}.class.name;
                Station{end+1,1} = R{i,k}.station.name;
                Rval(end+1) = RN(i,k);
            end
        end
    end
    Station = categorical(Station);
    JobClass = categorical(JobClass);    
    RespT = Rval(:); % we need to save first in a variable named like the column
    RT = Table(Station,JobClass,RespT);
    AvgTable = Table(Station,JobClass,RespT);
else
    Rval = zeros(M,K);
    JobClass = cell(K*M,1);
    Station = cell(K*M,1);
    for i=1:M
        for k=1:K
            JobClass{(i-1)*K+k} = R{i,k}.class.name;
            Station{(i-1)*K+k} = R{i,k}.station.name;
            Rval((i-1)*K+k) = RN(i,k);
        end
    end
    Station = categorical(Station);
    JobClass = categorical(JobClass);    
    RespT = Rval(:); % we need to save first in a variable named like the column
    RT = Table(Station,JobClass,RespT);
    AvgTable = Table(Station,JobClass,RespT);
end
end
