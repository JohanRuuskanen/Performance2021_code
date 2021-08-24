function [AvgSysChainTable, CT,XT] = getAvgSysTable(self,R,T)
% [AVGSYSCHAINTABLE, CT,XT] = GETAVGSYSTABLE(SELF,R,T)

% Return table of average system metrics
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.

if nargin==1
    R = self.model.getAvgRespTHandles;
    T = self.model.getAvgTputHandles;
end
[SysRespT, SysTput] = getAvgSys(self,R,T);
SysRespT = SysRespT';
SysTput = SysTput';
ChainObj = self.model.getChains();
Chain = cellfun(@(c) c.name,ChainObj,'UniformOutput',false)';
JobClasses = cell(0,1);
for c=1:length(Chain)    
    JobClasses(c,1) = {categorical(ChainObj{c}.classnames)};
end
Chain = categorical(Chain);
CT = Table(Chain, JobClasses, SysRespT);
XT = Table(Chain, JobClasses, SysTput);
AvgSysChainTable = Table(Chain, JobClasses, SysRespT, SysTput);
end
