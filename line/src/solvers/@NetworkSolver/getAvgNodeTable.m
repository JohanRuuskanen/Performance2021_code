function [AvgTable,QT,UT,RT,TT,AT] = getAvgNodeTable(self,Q,U,R,T,keepDisabled)
% [AVGTABLE,QT,UT,RT,TT] = GETNODEAVGTABLE(SELF,Q,U,R,T,KEEPDISABLED)
% Return table of average node metrics
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
if ~exist('keepDisabled','var')
    keepDisabled = false;
end
qn = self.model.getStruct;
I = self.model.getNumberOfNodes();
K = self.model.getNumberOfClasses();
if nargin == 1
    [Q,U,R,T,A] = self.model.getAvgHandles();
end

[QN,UN,RN,TN,AN] = self.getAvgNode(Q,U,R,T,A);

if isempty(QN)
    AvgTable = Table();
    QT = Table();
    UT = Table();
    RT = Table();
    TT = Table();
    AT = Table();
elseif ~keepDisabled
    Qval = []; Uval = [];
    Rval = []; Tval = []; Aval=[];
    JobClass = {};
    Node = {};
    for i=1:I
        for k=1:K
            if any(sum([QN(i,k),UN(i,k),RN(i,k),TN(i,k),AN(i,k)])>0)
                JobClass{end+1,1} = qn.classnames{k};
                Node{end+1,1} = qn.nodenames{i};
                Qval(end+1) = QN(i,k);
                Uval(end+1) = UN(i,k);
                Rval(end+1) = RN(i,k);
                Tval(end+1) = TN(i,k);
                Aval(end+1) = AN(i,k);
            end
        end
    end
    Node = categorical(Node);
    JobClass = categorical(JobClass);    
    QLen = Qval(:); % we need to save first in a variable named like the column
    QT = Table(Node,JobClass,QLen);
    Util = Uval(:); % we need to save first in a variable named like the column
    UT = Table(Node,JobClass,Util);
    RespT = Rval(:); % we need to save first in a variable named like the column
    RT = Table(Node,JobClass,RespT);
    Tput = Tval(:); % we need to save first in a variable named like the column
    TT = Table(Node,JobClass,Tput);
    ArvR = Aval(:); % we need to save first in a variable named like the column
    AT = Table(Node,JobClass,ArvR);
    AvgTable = Table(Node,JobClass,QLen,Util,RespT,ArvR,Tput);
else
    Qval = zeros(I,K); Uval = zeros(I,K);
    Rval = zeros(I,K); Tval = zeros(I,K); Aval = zeros(I,K);
    JobClass = cell(K*I,1);
    Node = cell(K*I,1);
    for i=1:I
        for k=1:K
            JobClass{end+1,1} = qn.classnames{k};
            Node{end+1,1} = qn.nodenames{i};
            Qval((i-1)*K+k) = QN(i,k);
            Uval((i-1)*K+k) = UN(i,k);
            Rval((i-1)*K+k) = RN(i,k);
            Tval((i-1)*K+k) = TN(i,k);
            Aval((i-1)*K+k) = AN(i,k);
        end
    end
    Node = categorical(Node);
    JobClass = categorical(JobClass);    
    QLen = Qval(:); % we need to save first in a variable named like the column
    QT = Table(Node,JobClass,QLen);
    Util = Uval(:); % we need to save first in a variable named like the column
    UT = Table(Node,JobClass,Util);
    RespT = Rval(:); % we need to save first in a variable named like the column
    RT = Table(Node,JobClass,RespT);
    Tput = Tval(:); % we need to save first in a variable named like the column
    TT = Table(Node,JobClass,Tput);
    ArvR = Aval(:); % we need to save first in a variable named like the column
    AT = Table(Node,JobClass,ArvR);
    AvgTable = Table(Node,JobClass,QLen,Util,RespT,ArvR,Tput);
end
end
