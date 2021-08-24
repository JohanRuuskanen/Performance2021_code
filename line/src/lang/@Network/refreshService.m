function [rates,scv,mu,phi,phases,lt] = refreshService(self)
% [RATES,SCV, MU,PHI,PHASES] = REFRESHSERVICE()
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
[rates, scv] = self.refreshRates;
[~,mu,phi,phases] = self.refreshServicePhases;
[lt] = self.refreshLST;

self.refreshScheduling(); % SEPT and LEPT may be affected
end
